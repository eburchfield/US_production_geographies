# Run all models
source("BS_load.R")
source("BS_clean.R")

# consider centering only
# https://stats.stackexchange.com/questions/29781/when-conducting-multiple-regression-when-should-you-center-your-predictor-varia

# This script runs all models (specification in run_model() function and renders RMDs for each crop.
corn <- run_model("corn", save=T)
soy <- run_model("soy", save=T)
wwheat <- run_model("wwheat", save=T)
cotton <- run_model("cotton", save=T)
hay <- run_model("hay", save=T)
alfalfa <- run_model("alfalfa", save=T)
sorghum <- run_model("sorghum", save=T)

rmarkdown::render("BS_analysis_template.RMD",
                  params = list(crop = "corn",
                                new_title = "Corn"),
                  output_file = "BS_analysis_corn_fert.html")
rmarkdown::render("BS_analysis_template.RMD",
                  params = list(crop = "alfalfa",
                                new_title = "Alfalfa"),
                  output_file = "BS_analysis_alfalfa_ns.html")
rmarkdown::render("BS_analysis_template.RMD",
                  params = list(crop = "cotton",
                                new_title = "Cotton"),
                  output_file = "BS_analysis_cotton_ns.html")
rmarkdown::render("BS_analysis_template.RMD",
                  params = list(crop = "sorghum",
                                new_title = "Sorghum"),
                  output_file = "BS_analysis_sorghum_ns.html")
rmarkdown::render("BS_analysis_template.RMD",
                  params = list(crop = "soy",
                                new_title = "Soy"),
                  output_file = "BS_analysis_soy_ns.html")
rmarkdown::render("BS_analysis_template.RMD",
                  params = list(crop = "hay",
                                new_title = "Hay"),
                  output_file = "BS_analysis_hay_ns.html")
rmarkdown::render("BS_analysis_template.RMD",
                  params = list(crop = "wwheat",
                                new_title = "Winter wheat"),
                  output_file = "BS_analysis_wwheat_ns.html")

rmarkdown::render("BS_model_comparisons.Rmd")
rmarkdown::render("BS_model_comparisons_ns.Rmd")
