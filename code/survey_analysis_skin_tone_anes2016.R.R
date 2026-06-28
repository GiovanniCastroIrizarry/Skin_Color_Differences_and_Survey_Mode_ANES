# ============================================================
# ANES 2016: Survey Mode and Self-Reported Skin Tone
# ============================================================
#
# Description:
# This script provides a reproducible survey-data analysis pipeline
# using the 2016 American National Election Study (ANES). The analysis
# evaluates whether self-reported skin tone varies by survey mode,
# comparing live face-to-face interviews and self-administered internet
# interviews.
#
# The workflow reads raw data directly from GitHub, merges respondent-
# level and interviewer-level files, constructs analysis-ready variables,
# performs data-quality checks through descriptive summaries and
# covariate balance diagnostics, estimates survey-weighted regression
# models, and generates figures for interpreting the results.
#
# Workflow:
#   1. Load respondent and methodology files from GitHub
#   2. Merge data sources using respondent identifiers
#   3. Clean and recode survey, demographic, political, and interviewer
#      variables
#   4. Produce descriptive statistics by survey mode and racial/ethnic group
#   5. Conduct two-sample t-tests comparing survey modes
#   6. Estimate survey-weighted OLS models for the full sample and
#      racial/ethnic subsamples
#   7. Generate coefficient plots, predicted values, and weighted
#      distribution figures
#   8. Assess covariate balance across survey modes
#   9. Estimate face-to-face robustness models using interviewer
#      characteristics
#
# Reproducibility:
#   - Data are read directly from GitHub raw URLs
#   - Tables are printed to the R console in markdown or tibble format
#   - Figures are displayed in the R plotting window
#
# Required packages:
#   readr, dplyr, tidyr, purrr, ggplot2, scales, ggtext, broom,
#   marginaleffects, modelsummary
#
# Data source:
#   American National Election Study (ANES) 2016 Time Series Study
#
# Repository:
#   https://github.com/GiovanniCastroIrizarry/Skin_Color_Differences_and_Survey_Mode_ANES
#
# ============================================================
# ------------------------------------------------------------
# Packages
# ------------------------------------------------------------
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(scales)
library(ggtext)
library(broom)
library(marginaleffects)
library(modelsummary)

# ------------------------------------------------------------
# 1. Load the data from GitHub
# ------------------------------------------------------------
github_base <- "https://raw.githubusercontent.com/GiovanniCastroIrizarry/Skin_Color_Differences_and_Survey_Mode_ANES/refs/heads/main/data/"

anes_main <- read_csv(paste0(github_base, "anes2016_main.csv"))
anes_methodology <- read_csv(paste0(github_base, "anes2016_methodology.csv"))

# Merge the interviewer variables onto the main file by
# respondent ID. A left join keeps every main-file respondent;
# the interviewer variables are only populated for face-to-face
# respondents (they are missing for self-administered cases).
anes_merged <- anes_main %>%
  left_join(anes_methodology, by = "V160001_orig")

# ------------------------------------------------------------
# 2. Build the analysis variables
# ------------------------------------------------------------
anes <- anes_merged %>%
  mutate(
    # Survey mode: factor with face-to-face (FTF) as the reference
    survey_mode = factor(
      V160501,
      levels = c(1, 2),
      labels = c("FTF", "Internet")
    ),
    self_administered = case_when(
      V160501 == 1 ~ 0,
      V160501 == 2 ~ 1,
      TRUE ~ NA_real_
    ),
    
    # Self-reported skin tone, 1 (lightest) to 10 (darkest)
    skin_tone_self = ifelse(V162368 < 0, NA, V162368),
    
    # Controls
    age    = ifelse(V161267 < 0, NA, V161267),
    gender = ifelse(V161342 < 0, NA, V161342),
    female = case_when(
      gender == 2 ~ 1,
      gender == 1 ~ 0,
      TRUE ~ NA_real_
    ),
    education = ifelse(V161270 < 0 | V161270 >= 90, NA, V161270),
    income    = ifelse(V161361x < 0, NA, V161361x),
    
    # Republican = 1, all other valid party ID = 0, DK/refused = NA
    republican = case_when(
      V161155 == 2 ~ 1,
      V161155 %in% c(0, 1, 3, 4) ~ 0,
      TRUE ~ NA_real_
    ),
    
    # Foreign-born indicator
    foreign_born = case_when(
      V161316 %in% c(1, 7) ~ 0,
      V161316 %in% c(2, 3, 4) ~ 1,
      TRUE ~ NA_real_
    ),
    
    # Race/ethnicity
    race = case_when(
      V161310x == 1 ~ "White",
      V161310x == 2 ~ "Black",
      V161310x == 3 ~ "Asian",
      V161310x == 5 ~ "Latino",
      TRUE ~ NA_character_
    ),
    
    # Interviewer-rated respondent skin tone, 1 (lightest) to 10 (darkest).
    # The ANES methodology file used the same Massey-Martin scale; observed
    # values top out below 10 in this sample, but 10 remains the valid ceiling.
    iwr_skin_tone = ifelse(V168302 >= 1 & V168302 <= 10, V168302, NA),
    iwr_hispanic = case_when(
      V168311 == 1 ~ 1,
      V168311 == 2 ~ 0,
      TRUE ~ NA_real_
    ),
    iwr_white = case_when(
      V168312 == 1 ~ 1,
      V168312 == 2 ~ 0,
      TRUE ~ NA_real_
    ),
    iwr_black = case_when(
      V168313 == 1 ~ 1,
      V168313 == 2 ~ 0,
      TRUE ~ NA_real_
    ),
    
    # Survey weight under a clear name
    weight = V160101
  ) %>%
  filter(!is.na(race))

# ------------------------------------------------------------
# 3. Descriptive statistics: skin tone by mode and racial group
# ------------------------------------------------------------
descriptive_table <- anes %>%
  filter(!is.na(skin_tone_self), !is.na(survey_mode)) %>%
  mutate(group = race) %>%
  bind_rows(
    anes %>%
      filter(!is.na(skin_tone_self), !is.na(survey_mode)) %>%
      mutate(group = "Full Sample")
  ) %>%
  group_by(group, survey_mode) %>%
  summarise(
    n    = n(),
    mean = round(mean(skin_tone_self, na.rm = TRUE), 2),
    sd   = round(sd(skin_tone_self, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  mutate(
    group = factor(
      group,
      levels = c("Full Sample", "White", "Latino", "Asian", "Black")
    )
  ) %>%
  arrange(group, survey_mode)

print(descriptive_table)

# ------------------------------------------------------------
# 4. Two-sample t-tests of skin tone by survey mode
#    Difference = Self-Administered (Internet) - Live Interview (FTF)
# ------------------------------------------------------------
add_stars <- function(p) {
  case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE ~ ""
  )
}

mode_ttest <- function(data) {
  d  <- data %>% filter(!is.na(skin_tone_self), !is.na(survey_mode))
  tt <- t.test(skin_tone_self ~ survey_mode, data = d)
  tibble(
    diff    = unname(tt$estimate[2] - tt$estimate[1]),  # Internet - FTF
    p_value = tt$p.value
  )
}

groups <- c("Full Sample", "White", "Latino", "Asian", "Black")

ttest_results <- map_dfr(groups, function(g) {
  d <- if (g == "Full Sample") anes else filter(anes, race == g)
  mode_ttest(d) %>% mutate(group = g, .before = 1)
}) %>%
  mutate(
    stars      = add_stars(p_value),
    difference = paste0(formatC(diff, format = "f", digits = 3), stars)
  )

print(ttest_results)

# ------------------------------------------------------------
# 5. Formatted descriptive table (Table 1) 
# ------------------------------------------------------------
table1 <- descriptive_table %>%
  mutate(cell = paste0(mean, " (", sd, ")\nN = ", n)) %>%
  select(group, survey_mode, cell) %>%
  pivot_wider(names_from = survey_mode, values_from = cell) %>%
  mutate(group = as.character(group)) %>%
  rename(
    `Racial Group`      = group,
    `Live Interview`    = FTF,
    `Self-Administered` = Internet
  )

diffs <- ttest_results %>%
  transmute(`Racial Group` = group, Difference = difference)

table1_final <- left_join(table1, diffs, by = "Racial Group")

print(table1_final)
#View(table1_final)

# ------------------------------------------------------------
# 6. Weighted OLS models of self-reported skin tone
# ------------------------------------------------------------

m1 <- lm(skin_tone_self ~ survey_mode + age + gender + education + income + republican,
         data = anes, weights = weight)

m2 <- lm(skin_tone_self ~ survey_mode + age + gender + education + income + republican,
         data = filter(anes, race == "White"), weights = weight)

m3 <- lm(skin_tone_self ~ survey_mode + age + gender + education + income + republican,
         data = filter(anes, race == "Black"), weights = weight)

m4 <- lm(skin_tone_self ~ survey_mode + age + gender + education + income + republican,
         data = filter(anes, race == "Asian"), weights = weight)

m5 <- lm(skin_tone_self ~ survey_mode + age + gender + education + income + republican,
         data = filter(anes, race == "Latino"), weights = weight)

m6 <- lm(skin_tone_self ~ survey_mode + age + gender + education + income + republican + foreign_born,
         data = filter(anes, race == "Asian"), weights = weight)

m7 <- lm(skin_tone_self ~ survey_mode + age + gender + education + income + republican + foreign_born,
         data = filter(anes, race == "Latino"), weights = weight)

# ------------------------------------------------------------
# 7. Survey-mode coefficient across all models (quick table)
# ------------------------------------------------------------
mode_effects <- map_dfr(
  list(
    `All Respondents`     = m1,
    White                 = m2,
    Black                 = m3,
    Asian                 = m4,
    Latino                = m5,
    `Asian (+ nativity)`  = m6,
    `Latino (+ nativity)` = m7
  ),
  ~ tidy(.x, conf.int = TRUE) %>% filter(term == "survey_modeInternet"),
  .id = "model"
) %>%
  select(model, estimate, std.error, statistic, p.value, conf.low, conf.high)

print(mode_effects)

# ------------------------------------------------------------
# 8. Table of OLS models
# ------------------------------------------------------------
modelsummary(
  list(
    "Full Sample"        = m1,
    "White"              = m2,
    "Black"              = m3,
    "Asian"              = m4,
    "Asian (+ nativity)" = m6,
    "Latino"             = m5,
    "Latino (+ nativity)" = m7
  ),
  coef_rename = c(
    "survey_modeInternet" = "Self-Administered mode",
    "age"                 = "Age",
    "gender"              = "Gender",
    "education"           = "Education",
    "income"              = "Income",
    "republican"          = "Republican",
    "foreign_born"        = "Foreign born"
  ),
  coef_omit = "(Intercept)",
  stars   = c("*" = 0.05, "**" = 0.01, "***" = 0.001),
  gof_map = c("nobs", "r.squared"),
  notes = paste(
    "OLS with survey weights. Controls include age, gender, education, income,",
    "and party identification. Models for the Asian and Latino subsamples",
    "additionally include a nativity (foreign-born) control. Standard errors",
    "in parentheses."
  ),
  output = "markdown"
)

# ------------------------------------------------------------
# 9. Figure: survey-mode coefficient across models
# ------------------------------------------------------------
coef_df <- map_dfr(
  list(
    Full_Sample    = m1,
    White          = m2,
    Black          = m3,
    Asian          = m4,
    Latino         = m5,
    Asian_nativity = m6,
    Latino_nativity = m7
  ),
  ~ tidy(.x, conf.int = TRUE) %>% filter(term == "survey_modeInternet"),
  .id = "model"
) %>%
  mutate(
    model = recode(
      model,
      "Full_Sample"     = "All Respondents",
      "Asian_nativity"  = "Asian\n(+ nativity)",
      "Latino_nativity" = "Latino\n(+ nativity)"
    ),
    model = factor(
      model,
      levels = c(
        "Latino\n(+ nativity)", "Latino",
        "Asian\n(+ nativity)", "Asian",
        "Black", "White", "All Respondents"
      )
    )
  )

fig_mode_effects <- ggplot(coef_df, aes(x = estimate, y = model,
                                        color = model, shape = model)) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.2, linewidth = 1) +
  geom_point(size = 3.5, stroke = 1) +
  scale_shape_manual(values = c(7, 6, 5, 4, 3, 2, 1)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Effect of Self-Administered (vs. Live Interview Mode)", y = "") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.y     = element_text(size = 11),
    panel.border    = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

fig_mode_effects

# ------------------------------------------------------------
# 10. Figure: predicted skin tone by survey mode and racial group
# ------------------------------------------------------------
preds_all <- predictions(
  m1, newdata = datagrid(survey_mode = c("FTF", "Internet"))
) %>%
  as.data.frame() %>%
  select(survey_mode, estimate, conf.low, conf.high) %>%
  mutate(group = "All\nRespondents")

preds_groups <- map_dfr(
  list(
    "White\n(non-Hispanic)" = m2,
    "Black"                 = m3,
    "Asian"                 = m4,
    "Latino"                = m5
  ),
  ~ predictions(.x, newdata = datagrid(survey_mode = c("FTF", "Internet"))) %>%
    as.data.frame() %>%
    select(survey_mode, estimate, conf.low, conf.high),
  .id = "group"
)

preds_all_df <- bind_rows(preds_all, preds_groups) %>%
  mutate(
    group = factor(
      group,
      levels = c("All\nRespondents", "White\n(non-Hispanic)",
                 "Latino", "Asian", "Black")
    ),
    survey_mode = factor(
      survey_mode,
      levels = c("FTF", "Internet"),
      labels = c("Live\nInterview", "Self\nAdministered")
    )
  )

fig_predictions <- ggplot(preds_all_df, aes(x = survey_mode, y = estimate, group = 1)) +
  geom_point(size = 4, color = "#2C3E50") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                width = 0.12, linewidth = 0.8, color = "#2C3E50") +
  facet_wrap(~ group, scales = "fixed", nrow = 1) +
  labs(x = "Survey Mode",
       y = "Predicted Skin Tone (1=Lightest, 10=Darkest)") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text         = element_text(face = "bold", size = 11),
    strip.background    = element_rect(fill = "gray95", color = NA),
    panel.border        = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text.x         = element_text(color = "black"),
    axis.text.y         = element_text(color = "black")
  )

fig_predictions

# ------------------------------------------------------------
# 11. Figure: skin tone distribution by survey mode (full sample)
# ------------------------------------------------------------
skin_dist <- anes %>%
  filter(!is.na(skin_tone_self), !is.na(survey_mode)) %>%
  group_by(survey_mode, skin_tone_self) %>%
  summarise(n = sum(weight), .groups = "drop") %>%
  group_by(survey_mode) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(survey_mode = factor(survey_mode,
                              levels = c("FTF", "Internet"),
                              labels = c("Live Interview", "Self-Administered")))

fig_skin_dist <- ggplot(skin_dist, aes(x = skin_tone_self, y = prop,
                                       color = survey_mode, linetype = survey_mode)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = 1:10) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_color_manual(values = c("Live Interview" = "#2C3E50",
                                "Self-Administered" = "#E74C3C")) +
  scale_linetype_manual(values = c("Live Interview" = "solid",
                                   "Self-Administered" = "dashed")) +
  labs(x = "Skin Tone (1 = Lightest, 10 = Darkest)",
       y = "Proportion of Respondents", color = NULL, linetype = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text        = element_text(color = "black")
  )

fig_skin_dist

# ------------------------------------------------------------
# 12. Figure: skin tone distribution by survey mode and racial group
# ------------------------------------------------------------
skin_dist_full <- anes %>%
  filter(!is.na(skin_tone_self), !is.na(survey_mode)) %>%
  group_by(survey_mode, skin_tone_self) %>%
  summarise(n = sum(weight), .groups = "drop") %>%
  group_by(survey_mode) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(
    survey_mode = factor(survey_mode,
                         levels = c("FTF", "Internet"),
                         labels = c("Live Interview", "Self-Administered")),
    group = "All Respondents"
  )

skin_dist_groups <- anes %>%
  filter(!is.na(skin_tone_self), !is.na(survey_mode), !is.na(race)) %>%
  group_by(race, survey_mode, skin_tone_self) %>%
  summarise(n = sum(weight), .groups = "drop") %>%
  group_by(race, survey_mode) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(
    survey_mode = factor(survey_mode,
                         levels = c("FTF", "Internet"),
                         labels = c("Live Interview", "Self-Administered")),
    group = race
  ) %>%
  select(-race)

skin_dist_all <- bind_rows(skin_dist_full, skin_dist_groups) %>%
  mutate(group = factor(group,
                        levels = c("All Respondents", "White",
                                   "Latino", "Asian", "Black")))

# Unweighted N per group for the facet labels
n_groups <- anes %>%
  filter(!is.na(skin_tone_self), !is.na(race)) %>%
  count(race, name = "n")

n_full <- anes %>%
  filter(!is.na(skin_tone_self)) %>%
  summarise(n = n()) %>%
  mutate(race = "All Respondents")

n_labels <- bind_rows(n_full, n_groups) %>%
  mutate(
    group = factor(race,
                   levels = c("All Respondents", "White",
                              "Latino", "Asian", "Black")),
    label = paste0(race, " (*N* = ", formatC(n, format = "d", big.mark = ","), ")")
  )

skin_dist_all <- skin_dist_all %>%
  left_join(n_labels %>% select(group, label), by = "group") %>%
  mutate(label = factor(label, levels = n_labels$label[order(n_labels$group)]))

fig_skin_dist_race <- ggplot(skin_dist_all, aes(x = skin_tone_self, y = prop,
                                                color = survey_mode, linetype = survey_mode)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = 1:10) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_color_manual(values = c("Live Interview" = "#2C3E50",
                                "Self-Administered" = "#E74C3C")) +
  scale_linetype_manual(values = c("Live Interview" = "solid",
                                   "Self-Administered" = "dashed")) +
  facet_wrap(~ label, ncol = 1, scales = "fixed") +
  labs(x = "Skin Tone (1 = Lightest, 10 = Darkest)",
       y = "Proportion of Respondents", color = NULL, linetype = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position  = "right",
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text        = element_text(color = "black"),
    strip.text       = element_markdown(face = "plain", size = 11),
    strip.background  = element_rect(fill = "gray95", color = NA)
  )

fig_skin_dist_race

# ------------------------------------------------------------
# 13. Covariate balance across survey modes (with standardized mean differences)
# ------------------------------------------------------------
balance_data <- anes %>%
  filter(
    !is.na(skin_tone_self), !is.na(survey_mode), !is.na(age),
    !is.na(gender), !is.na(education), !is.na(income),
    !is.na(republican), !is.na(weight)
  ) %>%
  mutate(
    mode_label = factor(survey_mode,
                        levels = c("FTF", "Internet"),
                        labels = c("Live Interview", "Self-Administered"))
  )

balance_table <- balance_data %>%
  select(mode_label, age, female, education, income, republican) %>%
  pivot_longer(-mode_label, names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  summarise(
    mean_live  = mean(value[mode_label == "Live Interview"], na.rm = TRUE),
    sd_live    = sd(value[mode_label == "Live Interview"], na.rm = TRUE),
    mean_self  = mean(value[mode_label == "Self-Administered"], na.rm = TRUE),
    sd_self    = sd(value[mode_label == "Self-Administered"], na.rm = TRUE),
    difference = mean_self - mean_live,
    pooled_sd  = sqrt((sd_live^2 + sd_self^2) / 2),
    smd        = difference / pooled_sd,
    p_value    = t.test(value ~ mode_label)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    Variable = recode(variable,
                      "age" = "Age", "female" = "Female",
                      "education" = "Education", "income" = "Income",
                      "republican" = "Republican"),
    `Live Interview`    = paste0(round(mean_live, 2), " (", round(sd_live, 2), ")"),
    `Self-Administered` = paste0(round(mean_self, 2), " (", round(sd_self, 2), ")"),
    Difference          = round(difference, 3),
    `Std. Difference`   = round(smd, 3),
    `p-value`           = round(p_value, 3)
  ) %>%
  select(Variable, `Live Interview`, `Self-Administered`,
         Difference, `Std. Difference`, `p-value`)

print(balance_table)
#View(balance_table)
datasummary_df(
  balance_table,
  title = "Covariate Balance Across Survey Modes",
  notes = paste(
    "Means (SD) by survey mode for the pooled analytic sample. Difference =",
    "Self-Administered minus Live Interview. Std. Difference is the difference",
    "divided by the pooled standard deviation. p-values from two-sample t-tests."
  ),
  output = "markdown"
)

# ------------------------------------------------------------
# 14. Figure: all coefficients by racial group
# ------------------------------------------------------------
coefplot_df <- map_dfr(
  list(
    "All Respondents|Base model" = m1,
    "White|Base model"           = m2,
    "Latino|Base model"          = m5,
    "Latino|+ nativity"          = m7,
    "Asian|Base model"           = m4,
    "Asian|+ nativity"           = m6,
    "Black|Base model"           = m3
  ),
  ~ tidy(.x, conf.int = TRUE),
  .id = "model_id"
) %>%
  separate(model_id, into = c("group", "model"), sep = "\\|") %>%
  filter(term != "(Intercept)") %>%
  mutate(
    term = recode(
      term,
      "survey_modeInternet" = "Self-administered mode",
      "age"                 = "Age",
      "gender"              = "Gender",
      "education"           = "Education",
      "income"              = "Income",
      "republican"          = "Republican",
      "foreign_born"        = "Foreign born"
    ),
    term = factor(
      term,
      levels = rev(c("Self-administered mode", "Age", "Gender",
                     "Education", "Income", "Republican", "Foreign born"))
    ),
    group = factor(group,
                   levels = c("All Respondents", "White",
                              "Latino", "Asian", "Black")),
    model = factor(model, levels = c("Base model", "+ nativity"))
  )

fig_coefplot_all <- ggplot(coefplot_df, aes(x = estimate, y = term,
                                            shape = model, linetype = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.6) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.18, linewidth = 0.8,
                 position = position_dodge(width = 0.55)) +
  geom_point(size = 2.8, stroke = 1, position = position_dodge(width = 0.55)) +
  facet_wrap(~ group, nrow = 1, scales = "fixed") +
  labs(x = "OLS Coefficient Estimate", y = "",
       shape = "Model", linetype = "Model") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position    = "bottom",
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    strip.text         = element_text(face = "bold", size = 11),
    strip.background    = element_rect(fill = "gray95", color = NA),
    panel.border        = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text.x         = element_text(color = "black"),
    axis.text.y         = element_text(color = "black")
  )

fig_coefplot_all

# ------------------------------------------------------------
# 15. Interviewer-race variables and robustness checks
#
# These analyses use the ANES methodology file to test whether the
# survey-mode effect is better explained by interviewer presence itself
# or by the racial composition of interviewers. The methodology file is
# necessary because interviewer race/ethnicity and interviewer-rated
# respondent skin tone are not in the main respondent file.
# ------------------------------------------------------------

anes <- anes %>%
  mutate(
    ftf = ifelse(self_administered == 0, 1, 0),
    
    # Mutually exclusive interviewer race category.
    # Reference category for models: Internet / no interviewer.
    iwr_race_cat = case_when(
      self_administered == 1 ~ "Internet",
      ftf == 1 & iwr_hispanic == 1 ~ "Hispanic",
      ftf == 1 & iwr_hispanic == 0 & iwr_black == 1 ~ "Black",
      ftf == 1 & iwr_hispanic == 0 & iwr_black == 0 & iwr_white == 1 ~ "White",
      ftf == 1 & iwr_hispanic == 0 & iwr_black == 0 & iwr_white == 0 ~ "Other",
      TRUE ~ NA_character_
    ),
    iwr_race_cat = factor(
      iwr_race_cat,
      levels = c("Internet", "Hispanic", "Black", "White", "Other")
    )
  )

# Sanity checks
print(table(anes$iwr_race_cat, anes$self_administered, useNA = "ifany"))
print(table(anes$iwr_race_cat, anes$race, useNA = "ifany"))


# ------------------------------------------------------------
# 16. Main models replacing survey-mode dummy with interviewer race
#
# Reference category: Internet / no interviewer.
# Each interviewer-race coefficient estimates how much darker respondents
# report their skin tone when interviewed by that category of interviewer
# relative to the self-administered internet condition.
# ------------------------------------------------------------

iwr_m1 <- lm(
  skin_tone_self ~ iwr_race_cat + age + gender + education + income + republican,
  data = anes,
  weights = weight
)

iwr_m2 <- lm(
  skin_tone_self ~ iwr_race_cat + age + gender + education + income + republican,
  data = filter(anes, race == "White"),
  weights = weight
)

iwr_m3 <- lm(
  skin_tone_self ~ iwr_race_cat + age + gender + education + income + republican,
  data = filter(anes, race == "Black"),
  weights = weight
)

iwr_m4 <- lm(
  skin_tone_self ~ iwr_race_cat + age + gender + education + income + republican,
  data = filter(anes, race == "Asian"),
  weights = weight
)

iwr_m5 <- lm(
  skin_tone_self ~ iwr_race_cat + age + gender + education + income + republican,
  data = filter(anes, race == "Latino"),
  weights = weight
)

iwr_m6 <- lm(
  skin_tone_self ~ iwr_race_cat + age + gender + education + income + republican + foreign_born,
  data = filter(anes, race == "Asian"),
  weights = weight
)

iwr_m7 <- lm(
  skin_tone_self ~ iwr_race_cat + age + gender + education + income + republican + foreign_born,
  data = filter(anes, race == "Latino"),
  weights = weight
)

# Table: interviewer-race models
modelsummary(
  list(
    "Full Sample"         = iwr_m1,
    "White"               = iwr_m2,
    "Black"               = iwr_m3,
    "Asian"               = iwr_m4,
    "Asian (+ nativity)"  = iwr_m6,
    "Latino"              = iwr_m5,
    "Latino (+ nativity)" = iwr_m7
  ),
  coef_rename = c(
    "iwr_race_catHispanic" = "Hispanic interviewer",
    "iwr_race_catBlack"    = "Black interviewer",
    "iwr_race_catWhite"    = "White interviewer",
    "iwr_race_catOther"    = "Other interviewer",
    "age"                  = "Age",
    "gender"               = "Gender",
    "education"            = "Education",
    "income"               = "Income",
    "republican"           = "Republican",
    "foreign_born"         = "Foreign born"
  ),
  coef_omit = "(Intercept)",
  stars   = c("*" = 0.05, "**" = 0.01, "***" = 0.001),
  gof_map = c("nobs", "r.squared"),
  notes = paste(
    "OLS with survey weights. Reference category is the internet condition,",
    "where no interviewer was present. Interviewer-race categories are mutually",
    "exclusive, with Hispanic ethnicity taking precedence."
  ),
  output = "markdown"
)

# Joint tests: do interviewer-race categories jointly improve model fit?
# This is the R equivalent of testing all interviewer-race coefficients together.
joint_iwr_test <- function(model) {
  reduced_model <- update(model, . ~ . - iwr_race_cat)
  anova(reduced_model, model)
}

cat("\nJoint F-test: Full sample\n")
print(joint_iwr_test(iwr_m1))

cat("\nJoint F-test: White respondents\n")
print(joint_iwr_test(iwr_m2))

cat("\nJoint F-test: Black respondents\n")
print(joint_iwr_test(iwr_m3))

cat("\nJoint F-test: Asian respondents\n")
print(joint_iwr_test(iwr_m4))

cat("\nJoint F-test: Latino respondents\n")
print(joint_iwr_test(iwr_m5))

cat("\nJoint F-test: Asian respondents + nativity\n")
print(joint_iwr_test(iwr_m6))

cat("\nJoint F-test: Latino respondents + nativity\n")
print(joint_iwr_test(iwr_m7))


# ------------------------------------------------------------
# 17. FTF-only interviewer race and racial concordance analyses
#
# These models restrict the sample to face-to-face respondents only.
# They test whether self-reported skin tone varies by interviewer
# race/ethnicity or by respondent-interviewer racial concordance.
# ------------------------------------------------------------

ftf <- anes %>% filter(ftf == 1)

# Descriptive: mean self-reported skin tone by respondent race and interviewer race
ftf_iwr_descriptives <- ftf %>%
  filter(!is.na(skin_tone_self), !is.na(iwr_race_cat)) %>%
  group_by(race, iwr_race_cat) %>%
  summarise(
    n = n(),
    weighted_n = sum(weight, na.rm = TRUE),
    mean_skin_tone = weighted.mean(skin_tone_self, weight, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(mean_skin_tone = round(mean_skin_tone, 2))

print(ftf_iwr_descriptives)

# Racial concordance indicator.
# Asian respondents are left missing because the methodology file does not
# contain a separate Asian-interviewer flag.
anes <- anes %>%
  mutate(
    concordance = case_when(
      ftf == 1 & race == "White"  & iwr_white == 1 ~ 1,
      ftf == 1 & race == "White"  & iwr_white == 0 ~ 0,
      ftf == 1 & race == "Black"  & iwr_black == 1 ~ 1,
      ftf == 1 & race == "Black"  & iwr_black == 0 ~ 0,
      ftf == 1 & race == "Latino" & iwr_hispanic == 1 ~ 1,
      ftf == 1 & race == "Latino" & iwr_hispanic == 0 ~ 0,
      TRUE ~ NA_real_
    ),
    concordance = factor(
      concordance,
      levels = c(0, 1),
      labels = c("Discordant", "Concordant")
    )
  )

ftf <- anes %>% filter(ftf == 1)

print(table(ftf$race, ftf$concordance, useNA = "ifany"))

conc_white <- lm(
  skin_tone_self ~ concordance + age + gender + education + income + republican,
  data = filter(ftf, race == "White"),
  weights = weight
)

conc_black <- lm(
  skin_tone_self ~ concordance + age + gender + education + income + republican,
  data = filter(ftf, race == "Black"),
  weights = weight
)

conc_latino <- lm(
  skin_tone_self ~ concordance + age + gender + education + income + republican + foreign_born,
  data = filter(ftf, race == "Latino"),
  weights = weight
)

modelsummary(
  list(
    "FTF White"  = conc_white,
    "FTF Black"  = conc_black,
    "FTF Latino" = conc_latino
  ),
  coef_rename = c(
    "concordanceConcordant" = "Concordant interviewer",
    "age"                   = "Age",
    "gender"                = "Gender",
    "education"             = "Education",
    "income"                = "Income",
    "republican"            = "Republican",
    "foreign_born"          = "Foreign born"
  ),
  coef_omit = "(Intercept)",
  stars   = c("*" = 0.05, "**" = 0.01, "***" = 0.001),
  gof_map = c("nobs", "r.squared"),
  notes = paste(
    "OLS with survey weights, face-to-face respondents only.",
    "Asian respondents are excluded because the ANES methodology file",
    "does not provide a separate Asian-interviewer indicator."
  ),
  output = "markdown"
)


# ------------------------------------------------------------
# 18. Self-reported vs. interviewer-rated skin tone, FTF only
#
# This is the most direct test of the mechanism. Among face-to-face
# respondents, both self-reported skin tone and interviewer-rated skin
# tone are available. If self-reports are anchored to observable
# phenotype under interviewer presence, the two measures should be
# strongly associated.
# ------------------------------------------------------------

ftf_skin_compare <- ftf %>%
  filter(!is.na(skin_tone_self), !is.na(iwr_skin_tone))

cat("\nSample size for FTF self-report vs. interviewer-rated comparison:\n")
print(nrow(ftf_skin_compare))

# Means side-by-side within FTF by respondent race
self_iwr_means <- ftf_skin_compare %>%
  group_by(race) %>%
  summarise(
    n = n(),
    self_mean = weighted.mean(skin_tone_self, weight, na.rm = TRUE),
    iwr_mean  = weighted.mean(iwr_skin_tone, weight, na.rm = TRUE),
    self_sd   = sd(skin_tone_self, na.rm = TRUE),
    iwr_sd    = sd(iwr_skin_tone, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    self_mean = round(self_mean, 2),
    iwr_mean  = round(iwr_mean, 2),
    self_sd   = round(self_sd, 2),
    iwr_sd    = round(iwr_sd, 2)
  )

print(self_iwr_means)

# Correlations: self-reported vs. interviewer-rated skin tone
cat("\nCorrelation: self-report vs. interviewer-rated skin tone, FTF only\n")
print(cor.test(ftf_skin_compare$skin_tone_self, ftf_skin_compare$iwr_skin_tone))

cat("\nCorrelations within racial/ethnic groups\n")
cor_results <- map_dfr(c("White", "Black", "Latino", "Asian"), function(g) {
  d <- ftf_skin_compare %>% filter(race == g)
  
  if (nrow(d) < 3) {
    return(tibble(group = g, n = nrow(d), correlation = NA_real_, p_value = NA_real_))
  }
  
  ct <- cor.test(d$skin_tone_self, d$iwr_skin_tone)
  
  tibble(
    group = g,
    n = nrow(d),
    correlation = unname(ct$estimate),
    p_value = ct$p.value
  )
})

print(cor_results)

# Signed gap: self-report minus interviewer rating.
# Negative values mean respondents self-report lighter than the interviewer rating.
ftf_skin_compare <- ftf_skin_compare %>%
  mutate(self_minus_iwr = skin_tone_self - iwr_skin_tone)

gap_table <- ftf_skin_compare %>%
  group_by(race) %>%
  summarise(
    n = n(),
    mean_gap = weighted.mean(self_minus_iwr, weight, na.rm = TRUE),
    sd_gap = sd(self_minus_iwr, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    mean_gap = round(mean_gap, 3),
    sd_gap = round(sd_gap, 3)
  )

print(gap_table)

# Paired t-tests: self-reported vs. interviewer-rated skin tone
cat("\nPaired t-test: all FTF respondents\n")
print(t.test(ftf_skin_compare$skin_tone_self, ftf_skin_compare$iwr_skin_tone, paired = TRUE))

cat("\nPaired t-tests by racial/ethnic group\n")
paired_tests <- map_dfr(c("White", "Black", "Latino", "Asian"), function(g) {
  d <- ftf_skin_compare %>% filter(race == g)
  
  if (nrow(d) < 3) {
    return(tibble(group = g, n = nrow(d), mean_difference = NA_real_, p_value = NA_real_))
  }
  
  tt <- t.test(d$skin_tone_self, d$iwr_skin_tone, paired = TRUE)
  
  tibble(
    group = g,
    n = nrow(d),
    mean_difference = unname(mean(d$skin_tone_self - d$iwr_skin_tone, na.rm = TRUE)),
    p_value = tt$p.value
  )
})

print(paired_tests)

# Regression: self-report as a function of interviewer-rated skin tone
anchor_all <- lm(
  skin_tone_self ~ iwr_skin_tone + age + gender + education + income + republican,
  data = ftf_skin_compare,
  weights = weight
)

anchor_white <- lm(
  skin_tone_self ~ iwr_skin_tone + age + gender + education + income + republican,
  data = filter(ftf_skin_compare, race == "White"),
  weights = weight
)

anchor_black <- lm(
  skin_tone_self ~ iwr_skin_tone + age + gender + education + income + republican,
  data = filter(ftf_skin_compare, race == "Black"),
  weights = weight
)

anchor_latino <- lm(
  skin_tone_self ~ iwr_skin_tone + age + gender + education + income + republican + foreign_born,
  data = filter(ftf_skin_compare, race == "Latino"),
  weights = weight
)

anchor_asian <- lm(
  skin_tone_self ~ iwr_skin_tone + age + gender + education + income + republican,
  data = filter(ftf_skin_compare, race == "Asian"),
  weights = weight
)

modelsummary(
  list(
    "FTF All"    = anchor_all,
    "FTF White"  = anchor_white,
    "FTF Black"  = anchor_black,
    "FTF Latino" = anchor_latino,
    "FTF Asian"  = anchor_asian
  ),
  coef_rename = c(
    "iwr_skin_tone" = "Interviewer skin tone rating",
    "age"           = "Age",
    "gender"        = "Gender",
    "education"     = "Education",
    "income"        = "Income",
    "republican"    = "Republican",
    "foreign_born"  = "Foreign born"
  ),
  coef_omit = "(Intercept)",
  stars   = c("*" = 0.05, "**" = 0.01, "***" = 0.001),
  gof_map = c("nobs", "r.squared"),
  notes = paste(
    "OLS with survey weights, face-to-face respondents only.",
    "The dependent variable is self-reported skin tone.",
    "The key predictor is interviewer-rated respondent skin tone."
  ),
  output = "markdown"
)


# ------------------------------------------------------------
# 19. Three-way comparison:
#     Internet self-report vs. FTF self-report vs. FTF interviewer rating
#
# This figure compares mean skin tone across three sources:
#   1. Internet self-report
#   2. FTF self-report
#   3. FTF interviewer rating
# ------------------------------------------------------------

weighted_sd <- function(x, w) {
  valid <- !is.na(x) & !is.na(w)
  x <- x[valid]
  w <- w[valid]
  
  if (length(x) <= 1) return(NA_real_)
  
  wm <- weighted.mean(x, w)
  sqrt(sum(w * (x - wm)^2) / sum(w))
}

weighted_summary <- function(data, variable, source_label) {
  x <- data[[variable]]
  w <- data$weight
  
  n <- sum(!is.na(x) & !is.na(w))
  m <- weighted.mean(x, w, na.rm = TRUE)
  s <- weighted_sd(x, w)
  se <- s / sqrt(n)
  
  tibble(
    source = source_label,
    mean = m,
    lb = m - 1.96 * se,
    ub = m + 1.96 * se,
    n = n
  )
}

triple_comparison <- map_dfr(c("White", "Black", "Asian", "Latino"), function(g) {
  internet_self <- anes %>%
    filter(race == g, self_administered == 1) %>%
    weighted_summary("skin_tone_self", "Internet self-report")
  
  ftf_self <- anes %>%
    filter(race == g, ftf == 1, !is.na(iwr_skin_tone)) %>%
    weighted_summary("skin_tone_self", "FTF self-report")
  
  ftf_iwr <- anes %>%
    filter(race == g, ftf == 1) %>%
    weighted_summary("iwr_skin_tone", "FTF interviewer rating")
  
  bind_rows(internet_self, ftf_self, ftf_iwr) %>%
    mutate(race = g)
}) %>%
  mutate(
    race = factor(race, levels = c("White", "Black", "Asian", "Latino")),
    source = factor(
      source,
      levels = c("Internet self-report", "FTF self-report", "FTF interviewer rating")
    )
  )

print(triple_comparison)

fig_triple_comparison <- ggplot(
  triple_comparison,
  aes(x = mean, y = source, group = 1)
) +
  geom_errorbar(
  aes(xmin = lb, xmax = ub),
  width = 0.12,
  linewidth = 0.8,
  orientation = "y"
) +
  geom_point(size = 3) +
  facet_wrap(~ race, ncol = 1) +
  labs(
    x = "Mean Skin Tone (1 = Lightest, 10 = Darkest)",
    y = ""
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black")
  )

fig_triple_comparison

# ============================================================
# End of script
# ============================================================