# Survival analysis for Broad-tailed Hummingbird at Rocky Mountain National Park
# Including age at first capture 

# Edited from original code by Erin Zylstra
# Gaby Samaniego
# gabysamaniego@arizona.edu
# 2025-06-11

# Load packages
library(tidyverse)
library(RMark)

# Clear environment
rm(list = ls()) 

# Load banding data
dat.raw <- read.csv('output/capture-data/cleanded-capture-data-RMNP-full.csv')

# Identify birds banded in September
sept.birds <- dat.raw %>% 
  filter(band_status == 1 & month == 9) %>% 
  select(band)

# Prepare data set for survival analysis 
dat <- dat.raw %>%
  filter(!band %in% sept.birds$band,# exclude September birds
         !band_site %in% c('WB2','WB1', 'WPK1', 'NFPC', 'POLC', 'SHIP'),
         month != 9) %>% # exclude September recaptures
  select(band, band_status, year, sex, obssite, band_age, band_site) %>% 
  rename(age = band_age) %>% 
  distinct() %>% 
  arrange(band, band_status, year)

# ------------------ Capture histories with age at first capture ------------- #
# ------------------------ includes juveniles and adults --------------------- #

# Create capture histories
ch.age <- dat %>%
  group_by(band, year, sex, age) %>%  
  summarize(n.observation = length(year), .groups = 'keep') %>%
  mutate(observed = 1) %>% 
  pivot_wider(names_from = year, 
              values_from = observed, 
              id_cols = c(band, sex, age), 
              values_fill = 0) %>% 
  relocate(band, '2003','2004','2005','2006','2007','2008','2009','2010','2011', 
           '2012', sex, age) %>%
  unite(ch, c('2003','2004','2005','2006','2007','2008','2009','2010', '2011','2012'), 
        sep = '') %>% 
  mutate(age = if_else(age == 'AHY', 'adult', 'young')) %>% 
  as.data.frame() %>% 
  select(-band)

# Make several variables factors (specifying levels for clarity)
ch.age$sex <- factor(ch.age$sex, levels = c('F', 'M'))
ch.age$age <- factor(ch.age$age, levels = c('young', 'adult'))

# Checks
head(ch.age)
str(ch.age)

# ------------------------------ PRPEPARE COVARIATES ------------------------- #

# Creates function to z-standardize the covariates
z.stand <- function(x) {
  (x - mean(x)) / sd(x)
}

# -------------------------------- Trapping Effort --------------------------- #

# Load effort data
effort.raw <- read.csv('output/banding-effort-data/banding-effort-all-sites-RMNP.csv')

# Edit effort data and standardize it
effort.z <- effort.raw %>% 
  group_by(year) %>%
  summarize(total_days = sum(total_banding_days, na.rm = TRUE),
            total_trap_hours = sum(total_trap_hours, na.rm = TRUE),
            .groups = 'drop') %>% 
  rename(time = year,
         effort_days = total_days,
         effort_hours = total_trap_hours) %>% 
  mutate(effort_days_z = z.stand(effort_days),
         effort_hours_z = z.stand(effort_hours)) %>%
  select(time, effort_days_z, effort_hours_z) %>% 
  as.data.frame()

# ---------------------------- Environmental Covariates ---------------------- # 

# Load data
winter.mx <- read.csv('output/weather-data/covariates-output/winter-covar-mexico.csv')
summer.co <- read.csv('output/weather-data/covariates-output/summer-covar-colorado.csv')

# Prepare covariates for analysis 
winter <- winter.mx %>% 
  mutate(time = 2002:2011, .after = winter_period) %>%
  select(time, aver_min_temp, aver_precip) %>% 
  mutate(winter_min_temp_z = z.stand(aver_min_temp),
         winter_precip_z = z.stand(aver_precip),
         winter_min_temp = aver_min_temp,
         winter_aver_precip = aver_precip, .keep = 'unused')

summer <- summer.co %>% 
  select(year, aver_max_temp, aver_min_temp) %>% 
  rename(time = year) %>% 
  mutate(summer_min_temp_z = z.stand(aver_min_temp),
         summer_max_temp_z = z.stand(aver_max_temp),
         summer_aver_min_temp = aver_min_temp,
         summer_aver_max_temp = aver_max_temp, .keep = 'unused')

# ----------------- PROCESS CAPTURE HISTORIES FOR MARK ANALYSIS -------------- #

# Process capture histories
age.process <- process.data(data = ch.age,
                            model = 'CJS',
                            begin.time = 2003,
                            groups = c('sex', 'age'),
                            age.var = 2,  
                            # Indicates that the second variable in 'groups' (age) 
                            # should be used to assign initial ages and track aging 
                            # over time
                            initial.ages = c(0, 1)) 
                            # Assigns starting age values for each level of 'age':
                            # 0 for birds first captured as juveniles
                            # 1 for birds first captured as adults

# Create design matrix
age.ddl <- make.design.data(age.process)

# This is the key step!
# Create age classes for survival where juvenile = age 0 and adults = 1+
age.ddl <- add.design.data(data = age.process,
                           ddl = age.ddl,
                           parameter = 'Phi',
                           type = 'age',
                           bins = c(0,1,9), 
                           # Creates 2 age classes using internal Age variable:
                           # [0, 1) for juveniles (Age = 0)
                           # [1, 9) for adults (Age = 1 to 8)
                           # Age = 9 is excluded with right = FALSE. This is fine
                           # as we don't have birds that survived from 2003 to 2012
                           right = FALSE, 
                           # If FALSE (default), intervals are left-closed, right-open:
                           # [a, b) includes a, excludes b
                           # If TRUE, intervals would be right-closed: (a, b]
                           name = 'ageclass') 
                           # The new variable 'ageclass' will be added to the design data
                           # and can now be used in model formulas like:
                           # Phi ~ ageclass

# Add effort to ddl 
age.ddl$p <- merge_design.covariates(age.ddl$p, effort.z)

# Add winter covariates to ddl
age.ddl$Phi <- merge_design.covariates(
  age.ddl$Phi, winter)

# Add summer covaraites to ddl
age.ddl$Phi <- merge_design.covariates(
  age.ddl$Phi, summer)

# Create a couple of other variables to help with model construction

# This code creates a new grouping variable called 'sexadult' in the survival
# design data. It defines three distinct groups:
  # Juveniles (0, both sexes combined). This will be the intercept when interpreting coefficients
  # Adult females (1)
  # Adult males (2)
age.ddl$Phi$sexadult <- ifelse(age.ddl$Phi$ageclass == '[0,1)', 0,
                               ifelse(age.ddl$Phi$sex == 'F', 1, 2))

# This categorical variable lets us model survival using a three-level factor.
# The model: Phi ~ sexadult translates to:
  # One survival estimate for juveniles
  # One for adult females
  # One for adult males

# Change new variable to a class factor
age.ddl$Phi$sexadult <- factor(age.ddl$Phi$sexadult)

# Create indicator variable for adults were Juvenile = 0 and Adult = 1
# I think this will be helpful for plots later on?
age.ddl$Phi$adult <- ifelse(age.ddl$Phi$ageclass == '[0,1)', 0, 1)

# Change new variable to a class factor
age.ddl$Phi$adult <- factor(age.ddl$Phi$adult)

# Inspect the newly created RMark objects
str(age.ddl)

# Phi
head(age.ddl$Phi, 10)
tail(age.ddl$Phi, 10)
summary(age.ddl$Phi)

# p
head(age.ddl$p, 10)
summary(age.ddl$p)

# ---------------------------------- RUN MODELS ------------------------------ #

# Run models exploring effects of sex and age class on survival in a 'building
# up model complexity' strategy

# -------------------------------- Using ageclass ---------------------------- #

# Is survival different between juveniles and adults?

# Create function
base.age.sex.models.1 <- function()
{
  Phi.dot <- list(formula = ~1) # based model
  Phi.age <- list(formula = ~ageclass) 
  
  p.sexEffort <- list(formula = ~effort_hours_z)
  
  cml <- create.model.list('CJS') 
  results <- mark.wrapper(cml, 
                          data = age.process,
                          ddl = age.ddl,
                          output = FALSE,
                          adjust = FALSE)
  return(results)
}

# Run function and store the results in a marklist
base.age.sex.results.1 <- base.age.sex.models.1()
base.age.sex.results.1

# Model with lowest Delta AIC
# Phi(~ageclass)p(~effort_hours_z) 0.0

# Look at estimates and standard errors 
results.1 <- base.age.sex.results.1[[1]]
results.1$results$beta

# The probability of survival differs by age class were adults have a significantly
# higher probability of survival than juveniles. 
# The probability of recapture increases with banding effort.

# The data supports the use of age as a variable in the survival models based on
# SE 

# Remove mark files so they don't clog repo
invisible(file.remove(list.files(pattern = 'mark.*\\.(inp|out|res|vcv|tmp)$')))

# -------------------------------- Using sexadult ---------------------------- #

# Is survival different among juveniles (females and males), adult females
# and adult males?

# Create function
base.age.sex.models.2 <- function()
{
  Phi.dot <- list(formula = ~1) 
  Phi.age <- list(formula = ~ageclass) 
  Phi.adultSex <- list(formula = ~sexadult)
  
  p.sexEffort <- list(formula = ~effort_hours_z)
  
  cml <- create.model.list('CJS') 
  results <- mark.wrapper(cml, 
                          data = age.process,
                          ddl = age.ddl,
                          output = FALSE,
                          adjust = FALSE)
  return(results)
}

# Run function and store the results in a marklist
base.age.sex.results.2 <- base.age.sex.models.2()
base.age.sex.results.2

# Model with lowest Delta AIC
# Phi(~sexadult)p(~effort_hours_z) 0.0
# Followed by far by 
# Phi(~ageclass)p(~effort_hours_z) 151.86

# Look at estimates and standard errors 
results.2 <- base.age.sex.results.2[[1]]
results.2$results$beta

# Juveniles (intercept) have a significantly lower probability (-) of survival 
# than both adult females and adult males. Among adults, females have a 
# significantly higher probability of survival than males.

# Remove mark files so they don't clog repo
invisible(file.remove(list.files(pattern = 'mark.*\\.(inp|out|res|vcv|tmp)$')))

# ----------------------------- Including covariates ------------------------- #

# -------------------- Adding wintering grounds covariates ------------------- #

# Is survival of juveniles and adults affected differently by the average minimum
# temperatures in Mexico?

# Create function
base.age.sex.models.3 <- function()
{
  Phi.sexAdult <- list(formula = ~sexadult)
  Phi.sexAdultCovar <- list(formula = ~sexadult + winter_min_temp_z)
  
  p.sexEffort <- list(formula = ~effort_hours_z)
  
  cml <- create.model.list('CJS') 
  results <- mark.wrapper(cml, 
                          data = age.process,
                          ddl = age.ddl,
                          output = FALSE,
                          adjust = FALSE)
  return(results)
}

# Run function and store the results in a marklist
base.age.sex.results.3 <- base.age.sex.models.3()
base.age.sex.results.3

# Model with lowest Delta AIC
# Phi(~sexadult + winter_min_temp_z)p(~effort_hours_z) 0.0
# Followed by far by 
# Phi(~sexadult)p(~sex + effort_hours_z) 34.32

# Look at estimates and standard errors 
results.3 <- base.age.sex.results.3[[2]]
results.3$results$beta

# Warmer winters increases the probability of survival of all groups 
# (+, significant)

# Look at real estimates
results.3$results$real

# Estimates for juveniles are around 0.2, for adult females are around 0.5
# and for adult males are around  0.3

# Remove mark files so they don't clog repo
invisible(file.remove(list.files(pattern = 'mark.*\\.(inp|out|res|vcv|tmp)$')))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Are juveniles more or less affected by the average min temperature in Mexico 
# than adult females and adult males? 

# Create function
base.age.sex.models.4 <- function()
{
  Phi.sexAdult <- list(formula = ~sexadult)
  Phi.sexAdultPluCovar <- list(formula = ~sexadult + winter_min_temp_z)
  Phi.sexAdultxCovar <- list(formula = ~sexadult * winter_min_temp_z)
  
  p.sexEffort <- list(formula = ~effort_hours_z)
  
  cml <- create.model.list('CJS') 
  results <- mark.wrapper(cml, 
                          data = age.process,
                          ddl = age.ddl,
                          output = FALSE,
                          adjust = FALSE)
  return(results)
}

# Run function and store the results in a marklist
base.age.sex.results.4 <- base.age.sex.models.4()
base.age.sex.results.4

# Model with lowest Delta AIC
# Phi(~sexadult * winter_min_temp_z)p(~effort_hours_z) 0.0
# Followed by
# Phi(~sexadult + winter_min_temp_z)p(~effort_hours_z) 27.28

# Look at estimates and standard errors 
results.4 <- base.age.sex.results.4[[3]]
results.4$results$beta

# The increase in min temperature is associated with increase in survival in juveniles 
# (+, not significant) 
# Adult females show a weaker response to the increase in temperature than juveniles,
# and seems that this increase has a negative effect on survival (-, not significant)
# Adult males show stronger positive response to the increase in temperature 
# compared to females and juveniles increasing survival (+, significant)

# Remove mark files so they don't clog repo
invisible(file.remove(list.files(pattern = 'mark.*\\.(inp|out|res|vcv|tmp)$')))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Is survival of juveniles and adults affected differently by the average 
# precipitation in Mexico 

# Create function
base.age.sex.models.5 <- function()
{
  Phi.sexAdult <- list(formula = ~sexadult)
  Phi.sexAdultPluCovar <- list(formula = ~sexadult + winter_precip_z)
  Phi.sexAdultxCovar <- list(formula = ~sexadult * winter_precip_z)
  
  p.sexEffort <- list(formula = ~effort_hours_z)
  
  cml <- create.model.list('CJS') 
  results <- mark.wrapper(cml, 
                          data = age.process,
                          ddl = age.ddl,
                          output = FALSE,
                          adjust = FALSE)
  return(results)
}

# Run function and store the results in a marklist
base.age.sex.results.5 <- base.age.sex.models.5()
base.age.sex.results.5

# Model with lowest Delta AIC
# Phi(~sexadult + winter_precip_z)p(~effort_hours_z) 0.0
# Followed by
# Phi(~sexadult * winter_precip_z)p(~effort_hours_z) 3.40

# Look at estimates and standard errors 
results.5 <- base.age.sex.results.5[[2]]
results.5$results$beta

# The increase in precipitation is associated with decrease survival for all groups
# (-, significant)

# Remove mark files so they don't clog repo
invisible(file.remove(list.files(pattern = 'mark.*\\.(inp|out|res|vcv|tmp)$')))


# ---------------------- Adding  summer grounds covariates ------------------- #

# Is survival of juveniles and adults affected differently by the average max
# temperature in summers in Colorado?

# Create function
base.age.sex.models.6 <- function()
{
  Phi.sexAdult <- list(formula = ~sexadult)
  Phi.sexAdultPluCovar <- list(formula = ~sexadult + summer_max_temp_z)
  Phi.sexAdultxCovar <- list(formula = ~sexadult * summer_max_temp_z)
  
  p.sexEffort <- list(formula = ~effort_hours_z)
  
  cml <- create.model.list('CJS') 
  results <- mark.wrapper(cml, 
                          data = age.process,
                          ddl = age.ddl,
                          output = FALSE,
                          adjust = FALSE)
  return(results)
}

# Run function and store the results in a marklist
base.age.sex.results.6 <- base.age.sex.models.6()
base.age.sex.results.6

# Model with lowest Delta AIC
# Phi(~sexadult + summer_max_temp_z)p(~effort_hours_z) 0.0
# Followed by 
# Phi(~sexadult * summer_max_temp_z)p(~effort_hours_z) 1.92

# Look at estimates and standard errors 
results.6 <- base.age.sex.results.6[[2]]
results.6$results$beta

# The increase in temperature is associated with an increase in the 
# probability of survival across all groups (+, significant). 

# When exploring the estimates for the model with the interaction, none were
# statistically significant

# Remove mark files so they don't clog repo
invisible(file.remove(list.files(pattern = 'mark.*\\.(inp|out|res|vcv|tmp)$')))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Is survival of juveniles and adults affected differently by the min temp in 
# summers in Colorado?

# Create function
base.age.sex.models.7 <- function()
{
  Phi.sexAdult <- list(formula = ~sexadult)
  Phi.sexAdultPluCovar <- list(formula = ~sexadult + summer_min_temp_z)
  Phi.sexAdultxCovar <- list(formula = ~sexadult * summer_min_temp_z)
  
  p.sexEffort <- list(formula = ~effort_hours_z)
  
  cml <- create.model.list('CJS') 
  results <- mark.wrapper(cml, 
                          data = age.process,
                          ddl = age.ddl,
                          output = FALSE,
                          adjust = FALSE)
  return(results)
}

# Run function and store the results in a marklist
base.age.sex.results.7 <- base.age.sex.models.7()
base.age.sex.results.7

# Model with lowest Delta AIC
# Phi(~sexadult + summer_min_temp_z)p(~effort_hours_z) 0.0
# Followed by 
# Phi(~sexadult * summer_min_temp_z)p(~effort_hours_z) 2.50

# Look at estimates and standard errors 
results.7 <- base.age.sex.results.7[[2]]
results.7$results$beta

# The increase in summer min temp (warmer nights) is associated with a decrease 
# in the probability of survival for all groups (-, significant). 

# When exploring the estimates for the model with the interaction, none were
# statistically significant

# Remove mark files so they don't clog repo
invisible(file.remove(list.files(pattern = 'mark.*\\.(inp|out|res|vcv|tmp)$')))


# ------------------------------- Run full model ----------------------------- #

# Create function
base.age.sex.full <- function() 
{
  Phi.full <- list(formula = ~sexadult * winter_min_temp_z 
                                       + winter_precip_z
                                       + summer_max_temp_z 
                                       + summer_min_temp_z)
  
  p.sexEffort <- list(formula = ~effort_hours_z)
  
  cml <- create.model.list('CJS') 
  results <- mark.wrapper(cml,
                          data = age.process,
                          ddl = age.ddl,
                          output = FALSE,
                          adjust = FALSE)
  return(results)
}

# Run function and store the results in a marklist
base.age.sex.full.results <- base.age.sex.full()

# Look at estimates
results.10 <- base.age.sex.full.results[[1]]
results.10$results$beta

# Phi:
# Juveniles (intercept) have lower probability of survival than adult males and 
# females (-, significant)
# Adult females have higher probability of survival than juveniles and males
# (+, significant)
# Adult males have higher probability of survival than juveniles but lower than 
# females (+, significant)

# The increase of precipitation in the wintering grounds increases the probability
# of survival for all groups, although the effect is small (+, barely not significant)
# The increase in max temperature (day temperature) in the summer grounds increases 
# the probability of survival for all groups, although the effect is small (+, significant)
# The increase in min temperature (night temperature) in the summer grounds 
# decreases the probability of survival for all groups (-, significant)

# Interaction:
# The increase in the min temp in the wintering grounds increases survival of 
# juveniles (+, barely not significant)
# For adult females, this increase has a small negative effect (-, not significant)
# The probability of survival of adult males increases as the min temperature in
# the wintering grounds increases (+, significant)

# p:
# More effort increases the probability of recapture (+, significant)

# Remove mark files so they don't clog repo
invisible(file.remove(list.files(pattern = 'mark.*\\.(inp|out|res|vcv|tmp)$')))


# --------------------------------- PLOT FINDINGS ---------------------------- #

# --------- Prepare data frames with real estimates for Phi and p ------------ #

# Extract real estimates and clean them up
real.ests <- results.10$results$real %>% 
  rownames_to_column(var = 'group') %>%
  mutate(parameter = ifelse(str_sub(group, 1, 3) == 'Phi', 'Phi', 'p'),
         group = ifelse(parameter == 'Phi', 
                        str_remove(group, 'Phi '), 
                        str_remove(group, 'p ')),
         year = as.numeric(str_sub(group, -4, -1)),
         a = str_sub(group, -7, -7),
         age = ifelse(a == '0', 'J', 'A'),
         sexMF = str_sub(group, 2, 2)) %>%
  select(parameter, year, age, sexMF, estimate, se, lcl, ucl)

# Extract real Phi estimates
real.Phi.ests <- real.ests %>%
  filter(parameter == 'Phi')

# Extract real p estimates
real.p.ests <- real.ests %>%
  filter(parameter == 'p') %>%
  select(-age)

# -------------------------------- Create plots ------------------------------ # 
# --------------------------------- Phi and p -------------------------------- #

# 1) Probability of survival of juveniles, adult females and adult males over time

# Prepare labels for plotting
real.Phi.ests <- real.Phi.ests %>%
  mutate(group = case_when(age == 'J' ~ 'Juvenile',
                           age == 'A' & sexMF == 'F' ~ 'Adult Female',
                           age == 'A' & sexMF == 'M' ~ 'Adult Male'))

# Plot survival probability
Phi.plot <- ggplot(real.Phi.ests, aes(x = as.numeric(year), 
                                      y = estimate,
                                      color = group,
                                      shape = group)) +
  geom_line(aes(group = group), size = 0.3) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0, linewidth = 0.3) +
  scale_color_manual(values = c('Juvenile' = 'gray10', 
                                'Adult Female' = 'gray26', 
                                'Adult Male' = 'gray46')) +
  scale_shape_manual(values = c('Juvenile' = 15,
                                'Adult Female' = 17,  
                                'Adult Male' = 19)) + 
  
  scale_y_continuous(limits = c(0,1), 
                     breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(breaks = 2003:2011) +
  labs(x = 'Year',
       y = 'Estimated annual survival probability\n(95% CI)',
       shape = 'Group',
       color = 'Group') +
  theme_classic() +
  theme(legend.position = 'right',
        plot.title = element_text(hjust = 0.5))
Phi.plot

# Save plot
ggsave(path = 'output/plots/New Plots Survival/',
       filename = 'survival.png',
       plot = Phi.plot,
       device = 'png',
       dpi = 300)

# 2) Recapture probability

# Plot recapture probability
p.plot <- ggplot(real.p.ests, aes(x = as.numeric(year), 
                                  y = estimate)) +
  geom_line(size = 0.3) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0, linewidth = 0.3) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(breaks = 2003:2012) +
  labs(y = 'Estimated recapture probability\n(95% CI)',
       x = 'Year') +
  theme_classic() +
  theme(legend.title = element_blank())
p.plot

# Save plot
ggsave(path = 'output/plots/New Plots Survival/',
       filename = 'recapture.png',
       plot = p.plot,
       device = 'png',
       dpi = 300)

# ---------------------------- Effect of covariates -------------------------- #

# Original base code from Erin Zylstra using BTLH data from Mount Lemmon
# Plots for NAEP poster 

# Code adapted and expanded using ChatGPT

# ------------------- First, prepare data needed for the plots --------------- #

# Merge covariables in a data frame
covars <- left_join(winter, summer, by = 'time')

# Extract Phi beta estimates
phi.betas <- results.10$results$beta %>%
  as.data.frame() %>%
  rownames_to_column('term') %>%
  filter(str_starts(term, 'Phi:')) %>%
  mutate(term = str_remove(term, 'Phi:'),
         term = str_replace(term, '\\(Intercept\\)', 'Intercept')) 

# Create a range of values for all covariates to plot

# Create function to build ranges for each covariate
build.range <- function(data, covar, n = 100) { # 100 seems standard? 
  values <- data[[covar]]  # Get the column named in 'covar' as a vector using [[ ]]
  seq(min(values), max(values), length.out = n)
}

# Create ranges 
winter.min.temp.range <- build.range(covars, 'winter_min_temp_z')
winter.precip.range <- build.range(covars, 'winter_precip_z')
summer.min.temp.range <- build.range(covars, 'summer_min_temp_z')
summer.max.temp.range <- build.range(covars, 'summer_max_temp_z')

# Create beta vectors from phi.betas
betas <- phi.betas$estimate

# Assign names to the beta vectors using the 'term' column in phi.betas
# This ensures that the coefficients match the correct columns in the 
# prediction data frames during matrix multiplication. 
names(betas) <- phi.betas$term

# To calculate confidence intervals for predictions by hand, we'll need the 
# variance-covariate matrix. Extract the values for Phi [1:9]
var.covar.matrix <- results.10$results$beta.vcv[1:9, 1:9]


# -------------------------------- Create plots ------------------------------ #

# 3) Effect of winter average min temperature on probability of survival on each group 
# (juveniles, adult females and adult males)

# Build prediction data frame
pred.df.winter.min.temp <- data.frame(
  Group = rep(c('Juvenile', 'Adult Female', 'Adult Male'), 
              each = length(winter.min.temp.range)),
  winter_min_temp_z = rep(winter.min.temp.range, times = 3)) # 3 for each group

# Add other covariates to predict data frame, holding them constant at 0 
# (standardized mean) 
pred.df.winter.min.temp$Intercept <- 1 
# Dummy variables for adult female and adult male
pred.df.winter.min.temp$sexadult1 <- ifelse(pred.df.winter.min.temp$Group == 'Adult Female', 1, 0)
pred.df.winter.min.temp$sexadult2 <- ifelse(pred.df.winter.min.temp$Group == 'Adult Male', 1, 0)
pred.df.winter.min.temp$summer_min_temp_z <- 0
pred.df.winter.min.temp$summer_max_temp_z <- 0
pred.df.winter.min.temp$winter_precip_z <- 0 

# Include interaction terms. Let the effect of temperature depend on sex
# Use `` (`sexadult2:winter_min_temp`) so the code works!
pred.df.winter.min.temp$`sexadult1:winter_min_temp_z` <- 
  pred.df.winter.min.temp$sexadult1 * pred.df.winter.min.temp$winter_min_temp_z
pred.df.winter.min.temp$`sexadult2:winter_min_temp_z` <- 
  pred.df.winter.min.temp$sexadult2 * pred.df.winter.min.temp$winter_min_temp_z 

# Calculate estimates of survival on the logit scale
# Selecting columns by names ensures that the order of variables in the 
# prediction data frame matches the order of the beta coefficients
estimate.winter.min.temp.logit <- as.matrix(pred.df.winter.min.temp[, names(betas)]) %*% 
  as.matrix(betas)

# Create design matrix
X.winter.min.temp <- as.matrix(pred.df.winter.min.temp[, names(betas)])

# Calculate bounds of confidence intervals on the logit scale
# A little matrix math, using the part of the var-covar matrix that applies to 
# survival parameters and not recapture parameters
std.errors.winter.min.temp <- sqrt(diag(X.winter.min.temp %*%
                                          var.covar.matrix %*%
                                          t(X.winter.min.temp)))
lcl.winter.min.temp.logit <- estimate.winter.min.temp.logit - 1.96 * std.errors.winter.min.temp
ucl.winter.min.temp.logit <- estimate.winter.min.temp.logit + 1.96 * std.errors.winter.min.temp

# Convert to probability scale
pred.df.winter.min.temp$estimate <- plogis(estimate.winter.min.temp.logit) 
pred.df.winter.min.temp$lcl <- plogis(lcl.winter.min.temp.logit)
pred.df.winter.min.temp$ucl <- plogis(ucl.winter.min.temp.logit)

# Back transform covariate to original scale
winter.min.temp.mean <- mean(covars$winter_min_temp)
winter.min.temp.sd <- sd(covars$winter_min_temp)
pred.df.winter.min.temp$winter_min_temp_c <- 
  pred.df.winter.min.temp$winter_min_temp_z * winter.min.temp.sd + winter.min.temp.mean  

# Plot
winter.min.temp.plot <- ggplot(pred.df.winter.min.temp, 
                               aes(x = winter_min_temp_c,
                                   y = estimate, 
                                   color = Group,
                                   fill = Group,
                                   linetype = Group)) +
  geom_line(size = 0.3) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.1, color = NA) +
  scale_x_continuous(breaks = seq(-2, 1, by = 0.5)) +
  labs(x = 'Average minimum temperature (°C) in wintering grounds',
       y = 'Estimated survival probability\n(95% CI)',
       color = 'Group',
       fill = 'Group',
       linetype = 'Group') +
  scale_color_manual(values = c('Juvenile' = 'gray10', 
                                'Adult Female' = 'gray26', 
                                'Adult Male' = 'gray46')) +
  scale_fill_manual(values = c('Juvenile' = 'gray10', 
                               'Adult Female' = 'gray26', 
                               'Adult Male' = 'gray46')) +
  scale_linetype_manual(values = c('Juvenile' = 'solid',
                                   'Adult Female' = 'longdash',
                                   'Adult Male' = 'dotted')) +
  theme_classic()
winter.min.temp.plot

# Save plot
ggsave(path = 'output/plots/New Plots Survival/',
       filename = 'winter min temp effect.png',
       plot = winter.min.temp.plot,
       device = 'png',
       dpi = 300)

# 4) Effect of summer min temp on probability of survival of all groups

# Build a prediction data frame
pred.df.summer.min.temp <- data.frame(
  Group = rep(c('Juvenile', 'Adult Female', 'Adult Male'), 
              each = length(summer.min.temp.range)),
  summer_min_temp_z = rep(summer.min.temp.range, times = 3)) 

# Add other covariates to predict data frame, holding them constant at 0 
# (standardized mean) 
pred.df.summer.min.temp$Intercept <- 1 
pred.df.summer.min.temp$sexadult1 <- ifelse(pred.df.summer.min.temp$Group == 'Adult Female', 1, 0)
pred.df.summer.min.temp$sexadult2 <- ifelse(pred.df.summer.min.temp$Group == 'Adult Male', 1, 0)
pred.df.summer.min.temp$winter_min_temp_z <- 0
pred.df.summer.min.temp$summer_max_temp_z <- 0
pred.df.summer.min.temp$winter_precip_z <- 0 
pred.df.summer.min.temp$`sexadult1:winter_min_temp_z` <- 0 # Use `` so the code works!
pred.df.summer.min.temp$`sexadult2:winter_min_temp_z` <- 0

# Calculate estimates of survival on the logit scale
# Selecting columns by names ensures that the order of variables in the 
# prediction data frame matches the order of the beta coefficients
estimate.summer.min.temp.logit <- as.matrix(pred.df.summer.min.temp[, names(betas)]) %*% 
  as.matrix(betas)

# Create design matrix
X.summer.min.temp <- as.matrix(pred.df.summer.min.temp[, names(betas)])

# Calculate bounds of confidence intervals on the logit scale
std.errors.summer.min.temp <- sqrt(diag(X.summer.min.temp %*%
                                          var.covar.matrix %*%
                                          t(X.summer.min.temp)))
lcl.summer.min.temp.logit <- estimate.summer.min.temp.logit - 1.96 * std.errors.summer.min.temp
ucl.summer.min.temp.logit <- estimate.summer.min.temp.logit + 1.96 * std.errors.summer.min.temp

# Convert to probability scale
pred.df.summer.min.temp$estimate <- plogis(estimate.summer.min.temp.logit) 
pred.df.summer.min.temp$lcl <- plogis(lcl.summer.min.temp.logit)
pred.df.summer.min.temp$ucl <- plogis(ucl.summer.min.temp.logit)

# Back transform covariate to original scale
summer.min.temp.mean <- mean(covars$summer_aver_min_temp)
summer.min.temp.sd <- sd(covars$summer_aver_min_temp)
pred.df.summer.min.temp$summer_min_temp_c <- 
  pred.df.summer.min.temp$summer_min_temp_z * summer.min.temp.sd + summer.min.temp.mean  

# Plot
summer.min.temp.plot <- ggplot(pred.df.summer.min.temp, 
                               aes(x = summer_min_temp_c,
                                   y = estimate, 
                                   color = Group,
                                   fill = Group,
                                   linetype = Group)) +
  geom_line(size = 0.3) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.1, color = NA) +
  labs(x = 'Average summer minimum temperature (°C)',
       y = 'Estimated survival probability\n(95% CI)',
       color = 'Group',
       fill = 'Group',
       linetype = 'Group') +
  scale_color_manual(values = c('Juvenile' = 'gray10', 
                                'Adult Female' = 'gray26', 
                                'Adult Male' = 'gray46')) +
  scale_fill_manual(values = c('Juvenile' = 'gray10', 
                               'Adult Female' = 'gray26', 
                               'Adult Male' = 'gray46')) +
  scale_linetype_manual(values = c('Juvenile' = 'solid',
                                   'Adult Female' = 'longdash',
                                   'Adult Male' = 'dotted')) +
  scale_x_continuous(breaks = seq(-12, -3.5, by = 1.5)) +
  theme_classic()
summer.min.temp.plot

# Save plot
ggsave(path = 'output/plots/New Plots Survival/', 
       filename = 'summer min temp effect.png',
       plot = summer.min.temp.plot,
       device = 'png',
       dpi = 300)

# 5) Effect of summer max temp on probability of survival of all groups

# Build a prediction data frame
pred.df.summer.max.temp <- data.frame(
  Group = rep(c('Juvenile', 'Adult Female', 'Adult Male'), 
              each = length(summer.max.temp.range)),
  summer_max_temp_z = rep(summer.max.temp.range, times = 3)) 

# Add other covariates to predict data frame, holding them constant at 0 
# (standardized mean) 
pred.df.summer.max.temp$Intercept <- 1 
pred.df.summer.max.temp$sexadult1 <- ifelse(pred.df.summer.max.temp$Group == 'Adult Female', 1, 0)
pred.df.summer.max.temp$sexadult2 <- ifelse(pred.df.summer.max.temp$Group == 'Adult Male', 1, 0)
pred.df.summer.max.temp$winter_min_temp_z <- 0
pred.df.summer.max.temp$summer_min_temp_z <- 0
pred.df.summer.max.temp$winter_precip_z <- 0 
pred.df.summer.max.temp$`sexadult1:winter_min_temp_z` <- 0 # Use `` so the code works!
pred.df.summer.max.temp$`sexadult2:winter_min_temp_z` <- 0

# Calculate estimates of survival on the logit scale
# Selecting columns by names ensures that the order of variables in the 
# prediction data frame matches the order of the beta coefficients
estimate.summer.max.temp.logit <- as.matrix(pred.df.summer.max.temp[, names(betas)]) %*% 
  as.matrix(betas)

# Create design matrix
X.summer.max.temp <- as.matrix(pred.df.summer.max.temp[, names(betas)])

# Calculate bounds of confidence intervals on the logit scale
std.errors.summer.max.temp <- sqrt(diag(X.summer.max.temp %*%
                                          var.covar.matrix %*%
                                          t(X.summer.max.temp)))
lcl.summer.max.temp.logit <- estimate.summer.max.temp.logit - 1.96 * std.errors.summer.max.temp
ucl.summer.max.temp.logit <- estimate.summer.max.temp.logit + 1.96 * std.errors.summer.max.temp

# Convert to probability scale
pred.df.summer.max.temp$estimate <- plogis(estimate.summer.max.temp.logit) 
pred.df.summer.max.temp$lcl <- plogis(lcl.summer.max.temp.logit)
pred.df.summer.max.temp$ucl <- plogis(ucl.summer.max.temp.logit)

# Back transform covariate to original scale
summer.max.temp.mean <- mean(covars$summer_aver_max_temp)
summer.max.temp.sd <- sd(covars$summer_aver_max_temp)
pred.df.summer.max.temp$summer_max_temp_c <- 
  pred.df.summer.max.temp$summer_max_temp_z * summer.max.temp.sd + summer.max.temp.mean  

# Plot
summer.max.temp.plot <- ggplot(pred.df.summer.max.temp, 
                               aes(x = summer_max_temp_c,
                                   y = estimate, 
                                   color = Group,
                                   fill = Group,
                                   linetype = Group)) +
  geom_line(size = 0.3) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.1, color = NA) +
  labs(x = 'Average summer maximum temperature (°C)',
       y = 'Estimated survival probability\n(95% CI)',
       color = 'Group',
       fill = 'Group',
       linetype = 'Group') +
  scale_color_manual(values = c('Juvenile' = 'gray10', 
                                'Adult Female' = 'gray26', 
                                'Adult Male' = 'gray46')) +
  scale_fill_manual(values = c('Juvenile' = 'gray10', 
                               'Adult Female' = 'gray26', 
                               'Adult Male' = 'gray46')) +
  scale_linetype_manual(values = c('Juvenile' = 'solid',
                                   'Adult Female' = 'longdash',
                                   'Adult Male' = 'dotted')) +
  scale_x_continuous(breaks = seq(27, 31, by = 0.5)) +
  theme_classic()
summer.max.temp.plot

# Save plot
ggsave(path = 'output/plots/New Plots Survival/', 
       filename = 'summer max temp effect.png',
       plot = summer.max.temp.plot,
       device = 'png',
       dpi = 300)

# 6) Effect of winter precip on probability of survival of all groups

# Build a prediction data frame
pred.df.winter.precip <- data.frame(
  Group = rep(c('Juvenile', 'Adult Female', 'Adult Male'), 
              each = length(winter.precip.range)),
  winter_precip_z = rep(winter.precip.range, times = 3))

# Add other covariates to predict data frame, holding them constant at 0 
# (standardized mean) 
pred.df.winter.precip$Intercept <- 1
pred.df.winter.precip$sexadult1 <- ifelse(pred.df.winter.precip$Group == 'Adult Female', 1, 0)
pred.df.winter.precip$sexadult2 <- ifelse(pred.df.winter.precip$Group == 'Adult Male', 1, 0)
pred.df.winter.precip$summer_min_temp_z <- 0
pred.df.winter.precip$winter_min_temp_z <- 0
pred.df.winter.precip$summer_max_temp_z <- 0
pred.df.winter.precip$`sexadult1:winter_min_temp_z` <- 0 
pred.df.winter.precip$`sexadult2:winter_min_temp_z` <- 0

# Calculate estimates of survival on the logit scale
estimate.winter.precip.logit <- as.matrix(pred.df.winter.precip[, names(betas)]) %*% 
  as.matrix(betas)

# Create design matrix
X.winter.precip <- as.matrix(pred.df.winter.precip[, names(betas)])

# Calculate bounds of confidence intervals on the logit scale
std.errors.winter.precip <- sqrt(diag(X.winter.precip %*%
                                        var.covar.matrix %*%
                                        t(X.winter.precip)))
lcl.winter.precip.logit <- estimate.winter.precip.logit - 1.96 * std.errors.winter.precip
ucl.winter.precip.logit <- estimate.winter.precip.logit + 1.96 * std.errors.winter.precip

# Convert to probability scale
pred.df.winter.precip$estimate <- plogis(estimate.winter.precip.logit) 
pred.df.winter.precip$lcl <- plogis(lcl.winter.precip.logit)
pred.df.winter.precip$ucl <- plogis(ucl.winter.precip.logit)

# Back transform covariate to original scale
winter.precip.mean <- mean(covars$winter_aver_precip)
winter.precip.sd <- sd(covars$winter_aver_precip)
pred.df.winter.precip$winter_aver_precip_c <- 
  pred.df.winter.precip$winter_precip_z * winter.precip.sd + winter.precip.mean  

# Plot
winter.precip.plot <- ggplot(pred.df.winter.precip, 
                             aes(x = winter_aver_precip_c,
                                 y = estimate, 
                                 color = Group,
                                 fill = Group,
                                 linetype = Group)) +
  geom_line(size = 0.3) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.1, color = NA) +
  labs(x = 'Winter precipitation (mm)',
       y = 'Estimated survival probability\n(95% CI)',
       color = 'Group',
       fill = 'Group',
       linetype = 'Group') +
  scale_color_manual(values = c('Juvenile' = 'gray10', 
                                'Adult Female' = 'gray26', 
                                'Adult Male' = 'gray46')) +
  scale_fill_manual(values = c('Juvenile' = 'gray10', 
                               'Adult Female' = 'gray26', 
                               'Adult Male' = 'gray46')) +
  scale_linetype_manual(values = c('Juvenile' = 'solid',
                                   'Adult Female' = 'longdash',
                                   'Adult Male' = 'dotted')) +
  theme_classic()
winter.precip.plot

# Save plot
ggsave(path = 'output/plots/New Plots Survival/',
       filename = 'winter precip effect.png',
       plot = winter.precip.plot,
       device = 'png',
       dpi = 300)
