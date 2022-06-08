## EEID 2022 Workshop Script
## by Noam Ross, MIT License (https://opensource.org/licenses/MIT)
## Find this script at https://bit.ly/eeid2022-wkshp

# Install these packages:
#install.packages(c(
#  "tidyverse",
#  "mgcv",
#  "sf",
#  "spdep",
#  "geofacet",
#  "lubridate",
#  "ggcorrplot",
#  "gratia"
#))
library(tidyverse)
library(mgcv)
library(spdep)
library(sf)
library(units)
library(geofacet)
library(lubridate)
library(broom)
library(ggcorrplot)
library(gratia)

## Download and import data.  CSVs are derived from the `tidycensus` package and
# https://health.data.ny.gov/Health/New-York-State-Statewide-COVID-19-Testing/xdss-u53e
nys_pop <- read_csv("https://dl.dropbox.com/s/a3v2ay7u0ma1nby/nys_pop.csv") %>%
  st_as_sf(wkt = "geometry", crs = 32618)

nys_cov <- read_csv("https://dl.dropbox.com/s/rvbv8xjfkjfauwp/nys_delta.csv")


#### Get to here? Great! If not, ask for help or buddy up with someone ####

# Let's look at our data. First, view and plot the raw data frames. The first
# `nys_cov`, is a subset of NY State COVID testing data during the delta surge
nys_cov
nys_pop
plot(nys_pop)

## Let's take a look at overall trends by aggregating the data and plotting
nys_all <- nys_cov %>%
  group_by(test_date) %>%
  summarize(new_positives = sum(new_positives),
            new_tests = sum(new_tests))

ggplot(nys_all, aes(x = test_date)) +
  geom_line(mapping = aes(y = new_positives))

ggplot(nys_all, aes(x = test_date)) +
  geom_line(mapping = aes(y = new_tests))

## Let's add population to our data frame and calculate some per-capita values
just_pop <- nys_pop %>%
  as_tibble() %>%
  select(-geometry)

nys_cov_pop <- nys_cov %>%
  left_join(just_pop, by = c("county")) %>%
  mutate(tests_per1000 = 1000*new_tests/population,
         pos_per1000 = 1000*new_positives/population,
         pop_dens = population / area_km2)

nys_cov_pop

## Now let's look at the variation of these values to think about how they
## affect our observations:

# Testing intensity:
nys_cov_pop %>% filter(test_date == as.Date("2020-12-01")) %>%
  ggplot(aes(x = tests_per1000)) +
  geom_histogram()

## What's that outlier?
nys_cov_pop %>% filter(test_date == as.Date("2020-12-01")) %>%
  arrange(desc(tests_per1000))

## Let's try the log scale
nys_cov_pop %>% filter(test_date == as.Date("2020-12-01")) %>%
  ggplot(aes(x = log10(tests_per1000))) +
  geom_histogram()

## Now let's look at other ways the population varies
nys_cov_pop %>% filter(test_date == as.Date("2020-12-01")) %>%
  ggplot(aes(x = log10(population))) +
  geom_histogram()

nys_cov_pop %>% filter(test_date == as.Date("2020-12-01")) %>%
  ggplot(aes(x = log10(pop_dens))) +
  geom_histogram()

## We need a good way to look at trends across the state, let's do something
## fun and create a geographic tiling.  The `geofacet` package does this.
## Check out https://cran.r-project.org/web/packages/geofacet/vignettes/geofacet.html

# Plot per-capita positive tests in counties across the state
nys_cov_pop %>%
  filter(test_date < as.Date("2020-11-01")) %>%
  ggplot(aes(x = test_date, y = pos_per1000)) +
  geom_line() +
  facet_geo(~county, grid = "us_ny_counties_grid1", scales = "free_y")

## So far we have learned that
## - Tests and cases have strong cyclic components
## - Areas vary hugely in population and density
## - Testing intensity varies hugely
## How can we make reasonable comparisons in time and space?
## Time to model!

# First, let's model our aggregate data, starting with a simple trend
# and then with a cyclic component

nys_all <- nys_all %>%
  mutate(wday = wday(test_date),
         day = as.integer(test_date - min(test_date)))
nys_all

mod_all0 <- gam(new_positives ~ s(day), # A simple single smooth
                data = nys_all,
                family = nb(),          # Use a negative binomial distribution
                method = "REML")  # Always pick this method, which will soon be default

summary(mod_all0)
plot(mod_all0)
gam.check(mod_all0)

preds <- predict(mod_all0, se.fit = TRUE) %>%
  as_tibble() %>%
  mutate(pred = exp(fit), lo = exp(fit - 2*se.fit), hi = exp(fit + 2*se.fit)) %>%
  bind_cols(nys_all)
preds
ggplot(preds, aes(x = day, y = new_positives)) +
  geom_line() +
  geom_line(mapping = aes(y = pred),
            col = 'red') +
  geom_ribbon(mapping = aes(ymin = lo, ymax = hi), fill="red", alpha = 0.25)

# OK, how to we deal with the cyclic affects? We add a "cc" - cyclic cubic - spline:
mod_all1 <- gam(new_positives ~ s(day, k=15) + s(wday, bs = "cc", k = 6),
                knots = list(wday=c(0.5, 7.5)),
                data = nys_all,
                family = nb(),
                method = "REML")

summary(mod_all1)
plot(mod_all1, pages = 1)

preds <- predict(mod_all1, se.fit = TRUE) %>%
  as_tibble() %>%
  mutate(pred = exp(fit), lo = exp(fit - 2*se.fit), hi = exp(fit + 2*se.fit)) %>%
  bind_cols(nys_all)
preds
ggplot(preds, aes(x = day, y = new_positives)) +
  geom_line() +
  geom_line(mapping = aes(y = pred),
            col = 'red') +
  geom_ribbon(mapping = aes(ymin = lo, ymax = hi), fill="red", alpha = 0.25)

## OK, Can we do better if we incorporate a measure of testing intensity?
mod_all2 <- gam(new_positives ~ s(day) + s(wday, bs = "cc", k = 6) +
                  s(new_tests),
                knots = list(wday=c(0.5, 7.5)),
                data = nys_all,
                family = nb(),
                method = "REML")
summary(mod_all2)
plot(mod_all2, pages = 1)

preds <- predict(mod_all2, se.fit = TRUE) %>%
  as_tibble() %>%
  mutate(pred = exp(fit), lo = exp(fit - 2*se.fit), hi = exp(fit + 2*se.fit)) %>%
  bind_cols(nys_all)
preds
ggplot(preds, aes(x = day, y = new_positives)) +
  geom_line() +
  geom_line(mapping = aes(y = pred),
            col = 'red') +
  geom_ribbon(mapping = aes(ymin = lo, ymax = hi), fill="red", alpha = 0.25)

# OK, cyclic effects and testing go a long way towards explaining data and
# maybe can help us make comparisons amongst our locations.  So let's use
# expand out .  We have a problem, though, spatial dependence!

# We can make a corrlation matrix of all the state measures:
cor_mat <- nys_cov %>%
  pivot_wider(id_cols = "test_date", names_from = "county", values_from = "new_positives") %>%
  select(-test_date) %>%
  cor()
ggcorrplot(cor_mat)

# In order to account for this, we're going to incorporate a "Markov Random Field",
# Where we explicitly account for correlation of adjacent counties.  Let's first
# For this, we have to create a neighborhood network, indicating which
# pairs of areas we'll model are adjacent.  Let's build this and visualize it:

nb <- poly2nb(nys_pop$geometry)
names(nb) <- nys_pop$county

centroids <- st_coordinates(st_centroid(nys_pop))
nb_net <-  nb2lines(nb = nb, coords = centroids, as_sf = TRUE)
st_crs(nb_net) <- st_crs(nys_pop)
ggplot() +
  geom_sf(data = nys_pop) +
  geom_sf(data = nb_net, col = "blue") +
  geom_point(data = as.data.frame(centroids), aes(x=X, y=Y), col = "blue")+
  theme_void()

# We now need to model a two-way interaction: the MRF in cases between counties,
# the change in cases over time, and their combination. `ti()` terms (tensor interactions)
# allow us to do this
nys_cov2 <- nys_cov %>%
  mutate(county = as.factor(county),  # gam() needs categories to be factors
         day = as.integer(test_date - min(test_date)),
         wday = wday(test_date)) %>%
  filter(test_date < as.Date("2020-12-01"))

mod_mrf1 <- gam(new_positives ~
                  ti(county, bs = "mrf", xt = list(nb = nb)) +
                  ti(day) +
                  ti(day, county, bs = c("ts", "mrf"), xt = list(nb = nb)),
                data = nys_cov2,
                family = nb(),
                method = "REML")
summary(mod_mrf1)

preds <- predict(mod_mrf1, se.fit = TRUE) %>%
  as_tibble() %>%
  mutate(pred = exp(fit), lo = exp(fit - 2*se.fit), hi = exp(fit + 2*se.fit)) %>%
  bind_cols(nys_cov2)
preds
ggplot(preds, aes(x = day, y = new_positives)) +
  geom_line() +
  geom_line(mapping = aes(y = pred),
            col = 'red') +
  geom_ribbon(mapping = aes(ymin = lo, ymax = hi), fill="red", alpha = 0.25) +
  facet_geo(~county, grid = "us_ny_counties_grid1", scales = "free_y")

#### Challenge time! ####

# 1) Incorporate the testing and day-of-week effects from our statewide model
#    into the county/MRF-level model.
# 2) View the model summary, the partial effects, and plot the geofaceted trends
# 3) Consider what variable makes for the most useful, comparable, and interpretable
#    effect of testing intensity - consider transforms, per-capita options, and try
#    these instead.
# 4) Bonus 1: Should testing or day-of-week effects vary by county? Why or why not?
#    Discuss and try incorporating this.
# 5) Bonus 2: How exactly does the MRF modify our predictions? Examine this by replacing
#    the MRF term with independent smooths for each county.  All three `ti()` terms
#    in the model can be replaced with: `s(day, by = county)` . Plot predictions and
#    compare how they behave compared to the MRF
#
#    (This can be slow to fit! Try replacing `gam` with `bam` ("big gam"),
#     use `method = "fREML`, ("fast REML"), and add the argument `discrete=TRUE``)


nys_cov3 <- nys_cov2 %>%
  left_join(as_tibble(nys_pop) %>% select(county, population),by = "county") %>%
  mutate(log_tests_per1000 = log10(new_tests/population),
         county = as.factor(county))

mod_mrf2 <- gam(new_positives ~
                  ti(county, bs = "mrf", xt = list(nb = nb)) +
                  ti(day) +
                  ti(day, county, bs = c("ts", "mrf"), xt = list(nb = nb)) +
                  s(wday, bs = "cc", k = 4) +
                  s(log_tests_per1000),
                knots = list(wday=c(0.5, 7.5)),
                data = nys_cov3,
                family = nb(),
                method = "REML")
summary(mod_mrf2)
plot(mod_mrf2, pages = 1)
preds <- predict(mod_mrf2, se.fit = TRUE) %>%
  as_tibble() %>%
  mutate(pred = exp(fit), lo = exp(fit - 2*se.fit), hi = exp(fit + 2*se.fit)) %>%
  bind_cols(nys_cov3)
preds
ggplot(preds, aes(x = day, y = new_positives)) +
  geom_line() +
  geom_line(mapping = aes(y = pred),
            col = 'red') +
  geom_ribbon(mapping = aes(ymin = lo, ymax = hi), fill="red", alpha = 0.25) +
  facet_geo(~county, grid = "us_ny_counties_grid1", scales = "free_y")

## Comparing non-MRF model:

mod_mrf3 <- bam(new_positives ~
                  s(day, by = county) +
                  s(wday, bs = "cc", k = 4) +
                  s(log_tests_per1000),
                knots = list(wday=c(0.5, 7.5)),
                data = nys_cov3,
                family = nb(),
                discrete = TRUE,
                method = "fREML")
summary(mod_mrf3)
plot(mod_mrf3)

preds_3 <- predict(mod_mrf3, se.fit = TRUE) %>%
  as_tibble() %>%
  mutate(pred = exp(fit), lo = exp(fit - 2*se.fit), hi = exp(fit + 2*se.fit)) %>%
  bind_cols(nys_cov3)
ggplot(preds_3, aes(x = day, y = new_positives)) +
  geom_line() +
  geom_line(mapping = aes(y = pred),
            col = 'red') +
  geom_line(data = preds, mapping = aes(y = pred), col = "blue") +
  facet_geo(~county, grid = "us_ny_counties_grid1", scales = "free_y")


## OK, we could keep going quite a bit, but now we need to figure out how to
## get stuff *out* of these models that is useful, namely:
## - Comparable values
## - Comporable growth rates

# We can make predictions for individual terms
terms <- predict(mod_mrf2, type = "terms", se.fit = TRUE)
terms
# But a simple way is just to replace terms with ones that make sense for us
# We can effectively remove the day-off week effects and normalize testing
# By setting these to constant values
nys_cov_norm <- nys_cov3 %>%
  mutate(wday = 4,
         log_tests_per1000 = quantile(log_tests_per1000, 0.9))

# Now make predictions with this normalized data frame and plot
preds <- predict(mod_mrf2, newdata = nys_cov_norm, se.fit = TRUE) %>%
  as_tibble() %>%
  mutate(pred = exp(fit), lo = exp(fit - 2*se.fit), hi = exp(fit + 2*se.fit)) %>%
  bind_cols(nys_cov_norm)
preds
ggplot(preds, aes(x = day, y = new_positives)) +
  geom_line() +
  geom_line(mapping = aes(y = pred),
            col = 'red') +
  geom_ribbon(mapping = aes(ymin = lo, ymax = hi), fill="red", alpha = 0.25) +
  facet_geo(~county, grid = "us_ny_counties_grid1", scales = "free_y")

# We are often interested in *rates* rather than values.  We can calculate
# derivatives by differences. For now I leave aside calculating the derivative
# variaces. The `gratia` simplifies this but does not yet handle models with MRFs,
# See the last example in ?predict.gam for manual calculation
eps = 1e-7 # a small finite difference
preds1 <- predict(mod_mrf2, newdata = nys_cov_norm)
preds2 <- predict(mod_mrf2, newdata = mutate(nys_cov_norm, day = day + eps))
preds <- preds %>%
  mutate(deriv = (preds2 - preds1) / eps)

ggplot(preds, aes(x = test_date, y = deriv)) +
  geom_line() +
  facet_geo(~county, grid = "us_ny_counties_grid1", scales = "fixed")

# OK, but this is the derivative in linear space. Since this is a log-scale
# model, it is equivalent to growth rate (r), in N = N0 * e^(r*t). We
# can use this to calculate case doubling time
preds <- preds %>%
  mutate(doubling_time = log(2) / deriv)

ggplot(preds, aes(x = test_date, y = doubling_time)) +
  geom_line() +
  facet_geo(~county, grid = "us_ny_counties_grid1", scales = "fixed")

## A (simplified) calculation of Rt can be derived from the growth rate
## and the mean generation (Tg) time of an infectious disease.
## Robust discussion in Walling & Lipsitch (2007, https://doi.org/10.1098%2Frspb.2006.3754)
# but one approximation R = 1 + r*Tg
##
Tg = 3 # An estimate of generation time for Delta COVID
preds <- preds %>%
  mutate(Rt = 1 + deriv*Tg)

ggplot(preds, aes(x = test_date, y = Rt)) +
  geom_line() +
  facet_geo(~county, grid = "us_ny_counties_grid1", scales = "fixed") +
  geom_hline(yintercept = 1, lty = 2, lwd = 0.25)

## Correlation among areas: Another thing that can be extracted from the MRF
## model is the degree of correlation among areas.  This is time-varying, so
## let's do it for a single point in time:

predmat <-  predict(mod_mrf2, type="lpmatrix", newdata = filter(nys_cov3, test_date == as.Date("2020-10-01")))
V <- vcov(mod_mrf2, unconditional = TRUE)
pred_vcov <- predmat %*% V %*% t(predmat)
pred_corr <- cov2cor(pred_vcov)

# Convert correlations to an sf network of lines with weights
nb_corr <- pred_corr * nb2mat(nb, style = "B", zero.policy = TRUE)
centroids <- st_coordinates(st_centroid(nys_pop))
nb_net <-  nb2lines(nb = nb, coords = centroids, as_sf = TRUE)
st_crs(nb_net) <- st_crs(nys_pop)
for (z in seq_len(nrow(nb_net))) {
  nb_net$wt[z] <- pred_corr[nb_net$i[z], nb_net$j[z]]
}

ggplot() +
  geom_sf(data = nys_pop) +
  geom_sf(data = nb_net, aes(col = wt, size = abs(wt))) +
  geom_point(data = as.data.frame(centroids), aes(x=X, y=Y)) +
  scale_color_distiller(type = "div", limits = c(-1,1)*max(abs(nb_net$wt)), name = "correlation") +
  scale_size_continuous(range = c(0.25, 2), name = "abs(correlation)") +
  theme_minimal()


