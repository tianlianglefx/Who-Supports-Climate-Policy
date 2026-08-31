# Who Supports Climate Policy?
# Reproducible analysis for ESS Round 8

# 1. Settings ---------------------------------------------------------------

SEED <- 1806
MAX_MSTOP <- 1500
B_CV <- 25
B_STAB <- 50
OUTPUT_DIR <- "outputs"

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

required_packages <- c(
  "readr", "dplyr", "ggplot2", "broom", "mboost", "stabs", "tibble"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install the following packages before running the analysis: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(readr)
library(dplyr)
library(ggplot2)
library(broom)
library(mboost)
library(stabs)
library(tibble)

# 2. Load data --------------------------------------------------------------

DATA_PATH <- "ESS8e02_3.csv"

if (!file.exists(DATA_PATH)) {
  stop(
    "ESS8e02_3.csv was not found in the working directory. ",
    "Download ESS Round 8 edition 2.3 from the ESS Data Portal and place the CSV here."
  )
}

ess_raw <- readr::read_csv(DATA_PATH, show_col_types = FALSE)

required_vars <- c(
  "cntry", "inctxff", "ccrdprs", "ccnthum", "ccgdbd", "lkredcc",
  "gvsrdcc", "agea", "gndr", "eisced", "hinctnta", "mnactic",
  "lrscale", "trstprl", "anweight"
)

missing_vars <- setdiff(required_vars, names(ess_raw))
if (length(missing_vars) > 0) {
  stop(
    "Required variables missing from the data: ",
    paste(missing_vars, collapse = ", ")
  )
}

# 3. Recode variables -------------------------------------------------------

valid_numeric <- function(x, valid_values) {
  ifelse(x %in% valid_values, as.numeric(x), NA_real_)
}

ess <- ess_raw %>%
  transmute(
    cntry = factor(cntry),

    inctxff_raw = valid_numeric(inctxff, 1:5),
    inctxff5 = factor(
      inctxff_raw,
      levels = 1:5,
      labels = c(
        "Strongly in favour",
        "Somewhat in favour",
        "Neither",
        "Somewhat against",
        "Strongly against"
      ),
      ordered = TRUE
    ),

    outcome_svo_num = case_when(
      inctxff %in% c(1, 2) ~ 1,
      inctxff %in% c(4, 5) ~ 0,
      TRUE ~ NA_real_
    ),
    outcome_svo_factor = factor(
      outcome_svo_num,
      levels = c(0, 1),
      labels = c("Opposition", "Support")
    ),

    outcome_svr_num = case_when(
      inctxff %in% c(1, 2) ~ 1,
      inctxff %in% c(3, 4, 5) ~ 0,
      TRUE ~ NA_real_
    ),
    outcome_svr_factor = factor(
      outcome_svr_num,
      levels = c(0, 1),
      labels = c("Rest", "Support")
    ),

    climate_skeptic = ccnthum == 55,
    structural_followup_na =
      ccrdprs == 66 | ccgdbd == 66 | lkredcc == 66 | gvsrdcc == 66,

    personal_responsibility = valid_numeric(ccrdprs, 0:10),
    climate_cause = case_when(
      ccnthum == 1  ~ "Entirely by natural processes",
      ccnthum == 2  ~ "Mainly by natural processes",
      ccnthum == 3  ~ "About equally by natural processes and human activity",
      ccnthum == 4  ~ "Mainly by human activity",
      ccnthum == 5  ~ "Entirely by human activity",
      ccnthum == 55 ~ "Does not think climate change is happening",
      TRUE ~ NA_character_
    ),
    climate_cause = factor(
      climate_cause,
      levels = c(
        "About equally by natural processes and human activity",
        "Entirely by natural processes",
        "Mainly by natural processes",
        "Mainly by human activity",
        "Entirely by human activity",
        "Does not think climate change is happening"
      )
    ),
    climate_harm = ifelse(
      ccgdbd %in% 0:10,
      10 - as.numeric(ccgdbd),
      NA_real_
    ),
    collective_efficacy = valid_numeric(lkredcc, 0:10),
    government_efficacy = valid_numeric(gvsrdcc, 0:10),

    age = ifelse(agea >= 15 & agea <= 100, as.numeric(agea), NA_real_),
    gender = case_when(
      gndr == 1 ~ "Male",
      gndr == 2 ~ "Female",
      TRUE ~ NA_character_
    ),
    gender = factor(gender, levels = c("Male", "Female")),
    education = ifelse(
      eisced %in% 1:7,
      as.character(eisced),
      NA_character_
    ),
    education = factor(education, levels = as.character(1:7)),
    income_decile = valid_numeric(hinctnta, 1:10),

    left_right = valid_numeric(lrscale, 0:10),
    trust_parliament = valid_numeric(trstprl, 0:10),

    main_activity = case_when(
      mnactic == 1 ~ "Paid work",
      mnactic == 2 ~ "Education",
      mnactic == 3 ~ "Unemployed, looking",
      mnactic == 4 ~ "Unemployed, not looking",
      mnactic == 5 ~ "Permanently sick/disabled",
      mnactic == 6 ~ "Retired",
      mnactic == 7 ~ "Community/military service",
      mnactic == 8 ~ "Housework/care",
      mnactic == 9 ~ "Other",
      TRUE ~ NA_character_
    ),
    main_activity = factor(
      main_activity,
      levels = c(
        "Paid work", "Education", "Unemployed, looking",
        "Unemployed, not looking", "Permanently sick/disabled", "Retired",
        "Community/military service", "Housework/care", "Other"
      )
    ),

    anweight = ifelse(
      is.finite(anweight) & anweight > 0,
      as.numeric(anweight),
      NA_real_
    )
  )

mean_anweight <- mean(ess$anweight, na.rm = TRUE)
if (!is.finite(mean_anweight) || mean_anweight <= 0) {
  stop("Could not construct model weights from anweight.")
}

ess <- ess %>%
  mutate(model_weight = anweight / mean_anweight)

# 4. Analysis samples -------------------------------------------------------

core_predictors <- c(
  "climate_cause", "age", "gender", "education", "income_decile",
  "left_right", "trust_parliament", "main_activity", "cntry"
)

full_predictors <- c(
  "personal_responsibility", "climate_cause", "climate_harm",
  "collective_efficacy", "government_efficacy", "age", "gender",
  "education", "income_decile", "left_right", "trust_parliament",
  "main_activity", "cntry"
)

make_model_sample <- function(data, outcome_num, outcome_factor, predictors) {
  required_for_model <- unique(c(
    outcome_num, outcome_factor, predictors, "model_weight"
  ))

  keep_for_diagnostics <- c(
    "anweight", "climate_skeptic", "structural_followup_na"
  )

  out <- data %>%
    select(any_of(unique(c(required_for_model, keep_for_diagnostics))))

  out[
    complete.cases(out[, required_for_model, drop = FALSE]),
    ,
    drop = FALSE
  ]
}

main_core <- make_model_sample(
  ess, "outcome_svo_num", "outcome_svo_factor", core_predictors
)
main_full <- make_model_sample(
  ess, "outcome_svo_num", "outcome_svo_factor", full_predictors
)
robust_core <- make_model_sample(
  ess, "outcome_svr_num", "outcome_svr_factor", core_predictors
)
robust_full <- make_model_sample(
  ess, "outcome_svr_num", "outcome_svr_factor", full_predictors
)

sample_flow <- tibble(
  step = c(
    "All respondents in pooled ESS8 sample",
    "Valid five-category fossil-fuel-tax outcome",
    "Main outcome eligible: support vs opposition",
    "Main core-model complete cases",
    "Main full-model complete cases",
    "Support-vs-rest outcome eligible",
    "Robustness core-model complete cases",
    "Robustness full-model complete cases"
  ),
  n = c(
    nrow(ess),
    sum(!is.na(ess$inctxff_raw)),
    sum(!is.na(ess$outcome_svo_num)),
    nrow(main_core),
    nrow(main_full),
    sum(!is.na(ess$outcome_svr_num)),
    nrow(robust_core),
    nrow(robust_full)
  )
)

write.csv(
  sample_flow,
  file.path(OUTPUT_DIR, "05_sample_flow.csv"),
  row.names = FALSE
)

make_structural_summary <- function(data, outcome_col, core_sample, full_sample, label) {
  eligible <- data[
    !is.na(data[[outcome_col]]) & !is.na(data$model_weight),
    ,
    drop = FALSE
  ]

  tibble(
    outcome_specification = label,
    eligible_n = nrow(eligible),
    eligible_structural_followup_na_n =
      sum(eligible$structural_followup_na, na.rm = TRUE),
    eligible_climate_skeptic_n =
      sum(eligible$climate_skeptic, na.rm = TRUE),
    core_complete_case_n = nrow(core_sample),
    core_climate_skeptic_n =
      sum(core_sample$climate_skeptic, na.rm = TRUE),
    full_complete_case_n = nrow(full_sample),
    full_climate_skeptic_n =
      sum(full_sample$climate_skeptic, na.rm = TRUE)
  )
}

structural_summary <- bind_rows(
  make_structural_summary(
    ess, "outcome_svo_num", main_core, main_full,
    "Main: support vs opposition"
  ),
  make_structural_summary(
    ess, "outcome_svr_num", robust_core, robust_full,
    "Robustness: support vs rest"
  )
)

write.csv(
  structural_summary,
  file.path(OUTPUT_DIR, "06_structural_missingness_summary.csv"),
  row.names = FALSE
)

# 5. Weighted descriptive statistics --------------------------------------

outcome_distribution <- ess %>%
  filter(!is.na(inctxff5), !is.na(anweight)) %>%
  group_by(inctxff5) %>%
  summarise(
    n = n(),
    weighted_n = sum(anweight),
    .groups = "drop"
  ) %>%
  mutate(weighted_percent = 100 * weighted_n / sum(weighted_n))

write.csv(
  outcome_distribution,
  file.path(OUTPUT_DIR, "07_weighted_outcome_distribution.csv"),
  row.names = FALSE
)

country_support <- ess %>%
  filter(!is.na(inctxff_raw), !is.na(anweight)) %>%
  group_by(cntry) %>%
  summarise(
    n_valid = n(),
    weighted_support_percent =
      100 * sum(anweight * as.numeric(inctxff_raw %in% c(1, 2))) / sum(anweight),
    .groups = "drop"
  ) %>%
  arrange(desc(weighted_support_percent))

write.csv(
  country_support,
  file.path(OUTPUT_DIR, "10_weighted_country_support.csv"),
  row.names = FALSE
)

p_country <- ggplot(
  country_support,
  aes(
    x = reorder(cntry, weighted_support_percent),
    y = weighted_support_percent
  )
) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Support for increasing fossil-fuel taxes by country",
    subtitle = "ESS analysis-weighted percentage choosing strongly/somewhat in favour",
    x = "Country",
    y = "Weighted percent supporting"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  file.path(OUTPUT_DIR, "11_weighted_country_support.png"),
  p_country,
  width = 8,
  height = 8,
  dpi = 300
)

# 6. Weighted logistic regression -----------------------------------------

formula_from_terms <- function(outcome, terms) {
  as.formula(paste(outcome, "~", paste(terms, collapse = " + ")))
}

main_core_formula <- formula_from_terms("outcome_svo_num", core_predictors)
main_full_formula <- formula_from_terms("outcome_svo_num", full_predictors)
robust_full_formula <- formula_from_terms("outcome_svr_num", full_predictors)

fit_weighted_glm <- function(formula, data) {
  glm(
    formula = formula,
    data = data,
    family = quasibinomial(link = "logit"),
    weights = model_weight
  )
}

baseline_main_core_same_sample <- fit_weighted_glm(
  main_core_formula,
  main_full
)
baseline_main_full <- fit_weighted_glm(
  main_full_formula,
  main_full
)
baseline_robust_full <- fit_weighted_glm(
  robust_full_formula,
  robust_full
)

tidy_or <- function(model) {
  broom::tidy(model) %>%
    mutate(
      OR = exp(estimate),
      conf.low = exp(estimate - 1.96 * std.error),
      conf.high = exp(estimate + 1.96 * std.error)
    ) %>%
    select(
      term, estimate, std.error, statistic, p.value,
      OR, conf.low, conf.high
    )
}

main_full_table <- tidy_or(baseline_main_full)
robust_full_table <- tidy_or(baseline_robust_full)

write.csv(
  main_full_table,
  file.path(OUTPUT_DIR, "14_main_full_weighted_OR.csv"),
  row.names = FALSE
)
write.csv(
  robust_full_table,
  file.path(OUTPUT_DIR, "16_support_vs_rest_full_weighted_OR.csv"),
  row.names = FALSE
)

weighted_fit_metrics <- function(model, data, outcome_col, weight_col = "model_weight") {
  y <- data[[outcome_col]]
  w <- data[[weight_col]]
  p <- predict(model, newdata = data, type = "response")

  eps <- 1e-12
  p <- pmin(pmax(p, eps), 1 - eps)

  weighted_log_loss <- -sum(
    w * (y * log(p) + (1 - y) * log(1 - p))
  ) / sum(w)

  p_null <- weighted.mean(y, w = w)
  p_null <- pmin(pmax(p_null, eps), 1 - eps)

  ll_model <- sum(w * (y * log(p) + (1 - y) * log(1 - p)))
  ll_null <- sum(w * (y * log(p_null) + (1 - y) * log(1 - p_null)))

  tibble(
    n = nrow(data),
    weighted_log_loss = weighted_log_loss,
    weighted_mcfadden_pseudo_R2 = 1 - ll_model / ll_null
  )
}

core_same_metrics <- weighted_fit_metrics(
  baseline_main_core_same_sample,
  main_full,
  "outcome_svo_num"
) %>%
  mutate(model = "Core model on full-model sample")

full_same_metrics <- weighted_fit_metrics(
  baseline_main_full,
  main_full,
  "outcome_svo_num"
) %>%
  mutate(model = "Full climate-attitude model")

core_full_fit_comparison <- bind_rows(
  core_same_metrics,
  full_same_metrics
) %>%
  select(model, n, weighted_log_loss, weighted_mcfadden_pseudo_R2)

write.csv(
  core_full_fit_comparison,
  file.path(OUTPUT_DIR, "17_main_core_vs_full_weighted_fit.csv"),
  row.names = FALSE
)

key_terms <- c(
  "personal_responsibility", "climate_harm", "collective_efficacy",
  "government_efficacy", "age", "income_decile", "left_right",
  "trust_parliament"
)

main_key_or <- main_full_table %>%
  filter(term %in% key_terms) %>%
  select(
    term,
    main_OR = OR,
    main_conf_low = conf.low,
    main_conf_high = conf.high,
    main_p_value = p.value
  )

robust_key_or <- robust_full_table %>%
  filter(term %in% key_terms) %>%
  select(
    term,
    support_vs_rest_OR = OR,
    support_vs_rest_conf_low = conf.low,
    support_vs_rest_conf_high = conf.high,
    support_vs_rest_p_value = p.value
  )

outcome_robustness_comparison <- full_join(
  main_key_or,
  robust_key_or,
  by = "term"
)

write.csv(
  outcome_robustness_comparison,
  file.path(OUTPUT_DIR, "18_outcome_robustness_weighted_OR_comparison.csv"),
  row.names = FALSE
)

# 7. Model-based boosting ---------------------------------------------------

continuous_boost_vars <- c(
  "personal_responsibility", "climate_harm", "collective_efficacy",
  "government_efficacy", "age", "income_decile", "left_right",
  "trust_parliament"
)

weighted_center <- function(x, w) {
  x - weighted.mean(x, w = w, na.rm = TRUE)
}

prepare_boost_data <- function(data) {
  out <- data
  for (v in continuous_boost_vars) {
    out[[paste0(v, "_c")]] <- weighted_center(
      out[[v]],
      out$model_weight
    )
  }
  out
}

main_boost_data <- prepare_boost_data(main_full)
robust_boost_data <- prepare_boost_data(robust_full)

boost_terms <- c(
  "bols(personal_responsibility_c, intercept = FALSE)",
  "bols(climate_harm_c, intercept = FALSE)",
  "bols(collective_efficacy_c, intercept = FALSE)",
  "bols(government_efficacy_c, intercept = FALSE)",
  "bols(left_right_c, intercept = FALSE)",
  "bols(trust_parliament_c, intercept = FALSE)",
  "bols(age_c, intercept = FALSE)",
  "bbs(age_c, center = TRUE, df = 1)",
  "bols(income_decile_c, intercept = FALSE)",
  "bbs(income_decile_c, center = TRUE, df = 1)",
  "bols(climate_cause, intercept = FALSE, df = 1)",
  "bols(gender, intercept = FALSE, df = 1)",
  "bols(education, intercept = FALSE, df = 1)",
  "bols(main_activity, intercept = FALSE, df = 1)",
  "brandom(cntry)"
)

make_boost_formula <- function(outcome_factor_col) {
  as.formula(
    paste(
      outcome_factor_col,
      "~",
      paste(boost_terms, collapse = " + ")
    )
  )
}

pretty_baselearner <- function(x) {
  case_when(
    grepl("personal_responsibility_c", x, fixed = TRUE) ~
      "Personal responsibility — linear",
    grepl("climate_harm_c", x, fixed = TRUE) ~
      "Climate harm — linear",
    grepl("collective_efficacy_c", x, fixed = TRUE) ~
      "Collective efficacy — linear",
    grepl("government_efficacy_c", x, fixed = TRUE) ~
      "Government efficacy — linear",
    grepl("left_right_c", x, fixed = TRUE) ~
      "Left-right orientation — linear",
    grepl("trust_parliament_c", x, fixed = TRUE) ~
      "Trust in parliament — linear",
    grepl("bbs(age_c", x, fixed = TRUE) ~
      "Age — smooth",
    grepl("bols(age_c", x, fixed = TRUE) ~
      "Age — linear",
    grepl("bbs(income_decile_c", x, fixed = TRUE) ~
      "Income decile — smooth",
    grepl("bols(income_decile_c", x, fixed = TRUE) ~
      "Income decile — linear",
    grepl("climate_cause", x, fixed = TRUE) ~
      "Climate cause — categorical",
    grepl("gender", x, fixed = TRUE) ~
      "Gender — categorical",
    grepl("education", x, fixed = TRUE) ~
      "Education — categorical",
    grepl("main_activity", x, fixed = TRUE) ~
      "Main activity — categorical",
    grepl("brandom(cntry)", x, fixed = TRUE) ~
      "Country — random intercept",
    TRUE ~ x
  )
}

run_weighted_boosting <- function(data, outcome_factor_col, tag, file_prefix) {
  boost_formula <- make_boost_formula(outcome_factor_col)

  set.seed(SEED)
  boost_model <- gamboost(
    formula = boost_formula,
    data = data,
    weights = data$model_weight,
    family = Binomial(link = "logit"),
    control = boost_control(
      mstop = MAX_MSTOP,
      nu = 0.1,
      trace = FALSE
    )
  )

  set.seed(SEED)
  subsampling_folds <- cv(
    model.weights(boost_model),
    type = "subsampling",
    B = B_CV,
    prob = 0.5
  )

  cv_risk <- cvrisk(
    boost_model,
    folds = subsampling_folds,
    grid = 1:MAX_MSTOP,
    papply = lapply
  )

  optimal_mstop <- mstop(cv_risk)
  if (optimal_mstop >= MAX_MSTOP) {
    warning("Optimal mstop reached the upper tuning boundary.")
  }

  capture.output(
    print(cv_risk),
    file = file.path(OUTPUT_DIR, paste0(file_prefix, "_cv_risk.txt"))
  )

  png(
    file.path(OUTPUT_DIR, paste0(file_prefix, "_cv_risk.png")),
    width = 1800,
    height = 1200,
    res = 180
  )
  plot(cv_risk)
  abline(v = optimal_mstop, lty = 2)
  dev.off()

  boost_tuned <- boost_model[optimal_mstop]

  selected_baselearners <- names(coef(boost_tuned))
  selected_table <- tibble(
    selected_baselearner = selected_baselearners,
    label = pretty_baselearner(selected_baselearners)
  )

  write.csv(
    selected_table,
    file.path(OUTPUT_DIR, paste0(file_prefix, "_selected_baselearners.csv")),
    row.names = FALSE
  )

  set.seed(SEED)
  stability_result <- stabsel(
    boost_tuned,
    q = 5,
    cutoff = 0.8,
    sampling.type = "SS",
    B = B_STAB,
    papply = lapply
  )

  capture.output(
    print(stability_result),
    file = file.path(OUTPUT_DIR, paste0(file_prefix, "_stability_selection.txt"))
  )

  stab_max <- NULL
  if (!is.null(stability_result$max)) {
    stab_max <- stability_result$max
  } else if (!is.null(stability_result$phat)) {
    phat <- stability_result$phat
    if (is.matrix(phat) || is.data.frame(phat)) {
      stab_max <- apply(phat, 1, max, na.rm = TRUE)
    } else {
      stab_max <- phat
    }
  }

  if (is.null(stab_max)) {
    stop("Could not extract stability-selection probabilities for ", tag, ".")
  }

  stability_table <- tibble(
    baselearner = names(stab_max),
    label = pretty_baselearner(names(stab_max)),
    selection_probability = as.numeric(stab_max)
  ) %>%
    arrange(desc(selection_probability))

  write.csv(
    stability_table,
    file.path(
      OUTPUT_DIR,
      paste0(file_prefix, "_stability_selection_probabilities.csv")
    ),
    row.names = FALSE
  )

  list(
    tag = tag,
    formula = boost_formula,
    tuned_model = boost_tuned,
    cv_risk = cv_risk,
    optimal_mstop = optimal_mstop,
    selected_table = selected_table,
    stability_table = stability_table
  )
}

main_boost_result <- run_weighted_boosting(
  data = main_boost_data,
  outcome_factor_col = "outcome_svo_factor",
  tag = "Main outcome: support vs opposition",
  file_prefix = "20_main_weighted_boost"
)

robust_boost_result <- run_weighted_boosting(
  data = robust_boost_data,
  outcome_factor_col = "outcome_svr_factor",
  tag = "Robustness outcome: support vs rest",
  file_prefix = "30_support_vs_rest_weighted_boost"
)

# 8. Comparison tables ------------------------------------------------------

boosting_tuning_summary <- tibble(
  outcome = c("Support vs opposition", "Support vs rest"),
  optimal_mstop = c(
    main_boost_result$optimal_mstop,
    robust_boost_result$optimal_mstop
  )
)

write.csv(
  boosting_tuning_summary,
  file.path(OUTPUT_DIR, "39_boosting_tuning_summary.csv"),
  row.names = FALSE
)

stability_comparison <- full_join(
  main_boost_result$stability_table %>%
    select(
      label,
      support_vs_opposition = selection_probability
    ),
  robust_boost_result$stability_table %>%
    select(
      label,
      support_vs_rest = selection_probability
    ),
  by = "label"
) %>%
  arrange(desc(support_vs_opposition), desc(support_vs_rest))

write.csv(
  stability_comparison,
  file.path(OUTPUT_DIR, "40_stability_selection_comparison.csv"),
  row.names = FALSE
)

key_stability_mapping <- tibble(
  term = c(
    "personal_responsibility", "climate_harm", "collective_efficacy",
    "government_efficacy", "age", "income_decile", "left_right",
    "trust_parliament"
  ),
  label = c(
    "Personal responsibility — linear",
    "Climate harm — linear",
    "Collective efficacy — linear",
    "Government efficacy — linear",
    "Age — linear",
    "Income decile — linear",
    "Left-right orientation — linear",
    "Trust in parliament — linear"
  )
)

logistic_vs_stability <- main_full_table %>%
  filter(term %in% key_stability_mapping$term) %>%
  select(
    term,
    weighted_logistic_OR = OR,
    conf.low,
    conf.high,
    p.value
  ) %>%
  left_join(key_stability_mapping, by = "term") %>%
  left_join(
    main_boost_result$stability_table %>%
      select(label, stability_selection_probability = selection_probability),
    by = "label"
  )

write.csv(
  logistic_vs_stability,
  file.path(OUTPUT_DIR, "41_key_predictors_logistic_vs_stability.csv"),
  row.names = FALSE
)

message("Analysis complete. Results saved to: ", normalizePath(OUTPUT_DIR))
