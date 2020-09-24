prep_data <- function(crop, scale = TRUE) {
  
  panel <- readRDS(paste0("./out/", crop, ".RDS"))
  panel <- panel %>% mutate(PERC_CROP = panel[,paste0(toupper(crop), "_SQKM")]/AG_SQKM)
  panel$PERC_AG <- panel$AG_SQKM/panel$TOTAL_SQKM
  
  clean_irr <- readRDS("C:/Users/eburchf/OneDrive - Emory University/WF/Diversity/out/clean_merged_irrigation.RDS") %>% 
    filter(YEAR %in% unique(panel$YEAR))
  
  # Variable selection based on expert opinion and collinearity
  panel <- panel %>%
    select(-c(UNITS, CROP,
              # remove land use metrics from other project
              SDI_CDL_BBOX_AG, SDI_CDL_AG, SDI_CDL_BBOX_ALL, SDI_CDL_ALL,
              SIDI_CDL_BBOX_AG, SIDI_CDL_AG, SIDI_CDL_BBOX_ALL, SIDI_CDL_ALL,
              RICH_CDL_AG, RICH_CDL_BBOX_AG, RICH_CDL_ALL, RICH_CDL_BBOX_ALL,
              RPR_CDL_BBOX_AG, RPR_CDL_AG, RPR_CDL_ALL, RPR_CDL_BBOX_ALL,
              # remove soil chars JC deemed nonessential (3)
              AWC_CLASS, ADD_PROP, T_GRAVEL, T_CEC_CLAY, T_CACO3, T_CASO4, REF_DEPTH,
              # remove irrigation term
              IRR_INT_SQKM, 
              # compute percent corn
              AG_SQKM,
              S_PH_H2O, # collinear with T_PH_H20
              T_BS, # collinear with T_PH_H2O
              T_TEB, # collinear with T_CEC_SOIL which was #1 from JC
              T_CLAY, # toss up between this and T_REF_BULK_DENSITY, went with latter b/c more broad
              LRR, # remove regions since this will be used in BS definition
              T_TEXTURE, T_ECE,
              DRAINAGE, # consistently low predictive power
              T_REF_BULK_DENSITY, T_SILT, T_SAND,
              TOTAL_SQKM, # KN EB 7/14
              PERC_IRR, # remove original, replace with clean version
              ELEVATION, SLOPE # highly collinear with region 7/14
    ))
  
  
  panel <- merge(panel, clean_irr, by = c("GEOID", "YEAR"))
  
  if (crop == "cotton") {
    panel$YIELD <- panel$YIELD/32
  }
  
  # drop rows with missing data (almost only yield)
  panel <- panel %>%
    na.omit()
  
  # drop regions with ten or fewer counties
  frr_cnt <- panel %>%
    group_by(FRR) %>%
    summarize(NG = length(unique(GEOID))) %>%
    filter(NG > 10)
  
  # drop counties more than one year available
  cty_cnt <- panel %>% 
    group_by(GEOID) %>% 
    count() %>% 
    filter(n > 1)
  
  # drop counties from panel
  panel <- panel %>% 
    filter(FRR %in% unique(frr_cnt$FRR)) %>%
    filter(GEOID %in% unique(cty_cnt$GEOID))
  
  # drop same counties from county shapefile
  county_sub <- county %>% 
    filter(GEOID %in% unique(panel$GEOID))
  county_sub <- county_sub %>% mutate(ID = 1:nrow(county_sub))
  
  county_df <- county_sub %>% select(-c(STATEFP, FRR))
  st_geometry(county_df) <- NULL
  panel <- merge(panel, county_df, by = "GEOID", all = T)

  if (scale == TRUE) {
      # scale predictors for comparison across models  
      scale <- panel %>% select(-c(ID, GEOID, YEAR, FRR)) 
      scale <- as.data.frame(scale(scale)) 
      out <- cbind.data.frame(panel %>% select(ID, GEOID, YEAR, FRR), scale)
  } else {
    
    out <- panel
  }
  outlist <- list(out, county_sub)
  return(outlist)
  
}

build_adj_matrix <- function(county_sub, crop, save=T) {

  temp <- poly2nb(county_sub, row.names = "ID", queen=T)
  h_adj <- nb2mat(temp, style="B", zero.policy = T) 
  h_adj <- as(h_adj, "dgTMatrix") # sparse style matrix conversion
  if (save == T) {
    saveRDS(h_adj, paste0("./out/hadj/hadj_", crop, ".RDS"))
  }
  return(h_adj)

}

run_model <- function(crop, save=T) {
  
  out <- prep_data(crop, scale=T)
  null <- out[[1]]
  county_sub <- out[[2]]
  hadj <- build_adj_matrix(county_sub, crop)

  null <- null %>% mutate_at(c("GDD", "TP"), ~inla.group(., n = 20, method = "quantile"))
  u <- sd(null$YIELD, na.rm = T)  
  
  apar <- 0.5 
  lmod <- lm(YIELD ~ GDD + TP  + T_CEC_SOIL + T_OC + T_PH_H2O + T_ESP + 
               PERC_IRR + 
               factor(YEAR), data=  null)
  bpar <- apar*var(residuals(lmod))
  lgprior <- list(prec = list(prior = "loggamma", param = c(apar, bpar))) 
  
  formula <- YIELD ~ 1 +  
    f(ID, model="bym2", 
      graph=hadj,
      scale.model=TRUE,
      constr = TRUE) +
    f(FRR, model = "iid", hyper = lgprior) +
    f(GDD, model = "rw2", 
      scale.model = T,
      hyper = list(theta = list(prior = "pc.prec", 
                                param = c(u, 0.01)))) +
    f(TP, model = "rw2", 
      scale.model = T,
      hyper = list(theta = list(prior = "pc.prec", 
                                param = c(u, 0.01)))) +
    PERC_IRR +
    T_CEC_SOIL + T_OC + T_PH_H2O + T_ESP + 
    factor(YEAR)
  
  print("Running model")
  mf <- inla(formula, data=null, family="gaussian",
             control.predictor = list(compute=T), 
             control.compute = list(dic=T, cpo=T))
  
  print("Diagnostics")
  df <- model_diagnostics(mf, null)
  
  out <- list(mf, df)
  
  if (save==T) {
    print("Saving")
    saveRDS(df, paste0("./out/models/df_", crop, ".RDS"))
    saveRDS(mf, paste0("./out/models/mf_", crop, ".RDS"))
  } 
  
  return(out)
}

model_diagnostics <- function(model_run, data) {
  
  if (is.character(model_run)) {
    r <- readRDS(model_run)
  } else (r <- model_run)
  
  fdf <- cbind(data, r$summary.fitted.values$mean)  # posterior distribution mean 
  colnames(fdf)[dim(fdf)[2]] <- c("FittedVals")
  
  # Deviance information criterion
  DIC <- r$dic$dic
  
  # Conditional predictive ordinate, LOOCV
  CPO <- sum(log(r$cpo$cpo), na.rm=T)
  
  # Probability Integral Transform
  PIT <- ggplot() +
    geom_histogram(aes(x = r$cpo$pit)) +
    ggtitle("PIT")
  
  # Posterior predictive distribution
  ppv <- c()
  for (i in 1:nrow(fdf)) {
    ppv[i] <- inla.pmarginal(q = fdf$YIELD[i], 
                             marginal = r$marginals.fitted.values[[i]])
  }
  
  # Linear indicates close fit between observed and predicted
  PPC_PT <- ggplot(fdf) +
    geom_point(aes(x = YIELD, y = fdf$FittedVals), alpha=0.2) +
    xlab("Observed") + ylab("Mean posterior predictive distribution") +
    theme_classic()
  
  # Reasonable fit should have few high P-values
  PPC_HIST <- ggplot() +
    geom_histogram(aes(x = ppv)) +
    ggtitle("Posterior predictive p-values")
  
  # MSE <- 1/(length(fdf$YIELD)*sum((fdf$YIELD - fdf$FittedVals)^2, na.rm=T))
  MSE <- MSE(fdf$FittedVals, fdf$YIELD)
  pred_res2<-(fdf$FittedVals[!is.na(fdf$YIELD)] - mean(fdf$YIELD, na.rm=T))^2
  obs_res2<-(fdf$YIELD[!is.na(fdf$YIELD)] - mean(fdf$YIELD, na.rm=T))^2
  R2<-sum(pred_res2, na.rm=T)/sum(obs_res2, na.rm=T)
  
  fit <- list(DIC, CPO, MSE, R2, PIT, PPC_PT, PPC_HIST)
  names(fit) <- c("DIC", "CPO", "MSE", "R2", "PIT", "PPC_PT", "PPC_HIST")
  return(fit)
  
}

total_area_effects <- function(model_run, data, level) {
  
  r <- model_run
  n <- length(unique(data[,level]))
  
  # extract area-specific residuals, random effects, Epsilon = u + v
  asr <- r$summary.random[[level]][1:n, c(1, 2, 3)]  # ID, mean, sd
  county_sp <- as(county_sub, "Spatial")
  
  if (level == "FRR") {
    asr <- r$summary.random[[level]][1:(n + 2), c(1, 2, 3)]  # ID, mean, sd
    colnames(asr) <- c("FRR", "mean", "SD")
  } 
  
  map_asr <- merge(county_sp, asr, by = level, all=T)
  map_asr <- merge(county, map_asr, by = 'GEOID', all = T)
  
  limit <- c(-1, 1)
  se <- ggplot() +
    geom_sf(data = map_asr, color = "transparent", size = 0.05, aes(fill = mean)) +
    geom_sf(data = state, color = "#A9A9A9", fill = NA, alpha = 0.1) +
    geom_sf(data = frr_shp, size = 1, fill = NA, alpha = 0.5) +
    project_theme +
    labs(fill = "") +
    scale_fill_gradient2(low = "#654321", mid = "#f5f5f5", high = "#000080", midpoint = 0, na.value = "#D3D3D3", 
                         limit = limit) +
    labs(fill = "",
         title = "",
         subtitle = "") +
    theme(text = element_text(size = 30), legend.text = element_text(size = 20)) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          rect = element_blank())
  
  
  return(se)
  
}

spatially_structured_effects <- function(model_run, data, level) {
  
  r <- model_run
  n <- length(unique(data[,level]))
  
  # extract spatially-structured residuals, u only
  asr <- r$summary.random[[level]][(n+1):(n*2), c(1, 2, 3)]  # ID, mean, sd
  asr$ID <- 1:n
  
  county_sp <- as(county_sub, "Spatial")
  
  map_asr <- sf::st_as_sf(merge(county_sp, asr, by = level))
  
  se <- ggplot(map_asr) +
    geom_sf(color = "transparent", size = 0.05, aes(fill = mean)) +
    geom_sf(data = state, color = "grey", fill = NA, alpha = 0.1) +
    theme_minimal() +
    geom_sf(data = frr_shp, size = 1, fill = "transparent") +
    project_theme +
    scale_fill_gradient2(low = "#d8b365", mid = "#f5f5f5", high = "#5ab4ac", midpoint = 0, na.value = "grey50")
  
  return(se)
  
}

spatially_unstructured_effects <- function(model_run, data, level) {
  
  r <- model_run
  n <- length(unique(data[,level]))
  
  # p. 187
  # "...the first n rows include area specific residuals (Epsilon = u+v), while the remaining present
  # information on the spatial structured residuals u_i only
  
  # extract spatially-structured residuals, u only
  asr_u <- r$summary.random[[level]][(n+1):(n*2), c(1, 2, 3)]  # ID, mean, sd
  asr_u$ID <- 1:n
  colnames(asr_u) <- c("ID", "MEAN_U", "SD_U")
  
  # extract area-specific residuals, random effects, Epsilon = u + v
  asr_uv <- r$summary.random[[level]][1:n, c(1, 2, 3)]  # ID, mean, sd
  colnames(asr_uv) <- c("ID", "MEAN_UV", "SD_UV")
  
  eff <- merge(asr_u, asr_uv, by = "ID")
  eff$V <- eff$MEAN_UV - eff$MEAN_U # this should isolate V only, or the unstructured residuals
  
  county_sp <- as(county_sub, "Spatial")
  
  map_asr <- sf::st_as_sf(merge(county_sp, eff, by = level))
  
  se <- ggplot(map_asr) +
    geom_sf(color = "transparent", size = 0.05, aes(fill =V)) +
    theme_minimal() +
    scale_fill_gradient2(low = "#d8b365", mid = "#f5f5f5", high = "#5ab4ac", midpoint = 0, na.value = "grey50") +
    geom_sf(data = state, color = "grey", fill = NA, alpha = 0.1) +
    geom_sf(data = frr_shp, size = 1, fill = "transparent") +
    project_theme
  
  out_list <- list(se, map_asr)
  
  return(out_list)
  
}

nonlinear_effect <- function(model_run, variable) {
  
  r <- model_run
  nle <- ggplot(data = r$summary.random[[variable]][,c(1,4:6)], 
                aes(x = ID, y = `0.5quant`)) +
    geom_line() +
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), alpha = 0.3) +
    theme_minimal() +
    theme(legend.position="none") +
    xlab(variable) +
    ylab("Effect on yield (bu/ac)")
  
  return(nle)
}

plot_coef <- function(model) {
  
  coef <- model$summary.fixed[,c(1, 3, 5)] # mean, 0.025, 0.985
  colnames(coef) <- c("MEAN", "LOW", "HIGH")
  coef$VARS <- rownames(coef)
  coef <- coef %>% filter(VARS != "(Intercept)")
  # drop year factors for now and intercept
  # coef <- coef[2:9,]
  
  ggplot(coef) +
    geom_errorbarh(height = 0, aes(xmin = LOW, xmax = HIGH, 
                                   y = reorder(VARS, desc(LOW)))) +
    project_theme +
    geom_point(aes(x = MEAN, y = reorder(VARS, desc(LOW))), size = 2) +
    geom_vline(xintercept = 0) +
    labs(x = "",
         y = "")
}

compute_diff <- function(model, data, level) {
  
  r <- model
  b0 <- r$summary.fixed[1,1]
  
  # extract area-specific residuals
  n <- length(unique(data[,"ID"]))
  asr <- r$summary.random[["ID"]][1:n, c(1, 2)]  # E_ij = u_0jk + v_0jk
  colnames(asr) <- c("ID", "COUNTY_MEAN")
  
  # extract area-specific residuals
  n2 <- length(unique(data[,level]))
  asr2 <- r$summary.random[[level]][ , c(1, 2)]  # FRR, mean
  asr2$mean <- asr2$mean # v_00k
  colnames(asr2) <- c("FRR", "REGION_MEAN")
  
  m <- merge(county_sub, asr, by = "ID", all=T)
  m <- merge(m, asr2, by = level, all=T)
  m$DIFF <- m$COUNTY_MEAN - m$REGION_MEAN
  m$Z <- NA
  
  lrr_ids <- unique(asr2[,level])
  for (i in 1:length(lrr_ids)) {
    id <- lrr_ids[[i]]
    df <- m %>% filter(FRR == id)
    df$Z <- (df$DIFF - mean(df$DIFF, na.rm=T))/sd(df$DIFF, na.rm=T)
    m$Z[m$FRR == id] <- df$Z
  }
  
  return(m)
  
}

compute_diff_v <- function(model, data, level) {
  
  r <- model
  b0 <- r$summary.fixed[1,1]
  n <- length(unique(data[,"ID"]))
  
  # extract spatially-structured residuals, u only
  asr_u <- r$summary.random[["ID"]][(n+1):(n*2), c(1, 2, 3)]  # ID, mean, sd
  asr_u$ID <- 1:n
  colnames(asr_u) <- c("ID", "MEAN_U", "SD_U")
  
  # extract area-specific residuals, random effects, Epsilon = u + v
  asr_uv <- r$summary.random[["ID"]][1:n, c(1, 2, 3)]  # ID, mean, sd
  colnames(asr_uv) <- c("ID", "MEAN_E", "SD_E")
  
  eff <- merge(asr_u, asr_uv, by = "ID")
  eff$MEAN_V <- eff$MEAN_E - eff$MEAN_U # this should isolate V only, or the unstructured residuals
  
  # extract area-specific residuals
  n2 <- length(unique(data[,level]))
  asr2 <- r$summary.random[[level]][ , c(1, 2)]  # FRR, mean
  asr2$mean <- asr2$mean # v_00k
  colnames(asr2) <- c("FRR", "REGION_MEAN")
  
  m <- merge(county_sub, eff, by = "ID", all=T)
  m <- merge(m, asr2, by = level, all=T)
  
  m$DIFF <- m$MEAN_V - m$REGION_MEAN # county - region, so high diff means county has higher yields, negative means region has higher yields
  m$Z <- NA
  
  lrr_ids <- unique(asr2[,level])
  for (i in 1:length(lrr_ids)) {
    id <- lrr_ids[[i]]
    df <- m %>% filter(FRR == id)
    df$Z <- (df$DIFF - mean(df$DIFF, na.rm=T))/sd(df$DIFF, na.rm=T)
    m$Z[m$FRR == id] <- df$Z
  }
  
  return(m)
  
}

make_shp <- function(mlist, crops, save=F) {
  
  cty <- county
  
  for (i in 1:length(mlist)) {
    
    df <- prep_data(crops[i])[[1]]
    county_sub <- prep_data(crops[i])[[2]]
    mc <- compute_diff(mlist[[i]], df, level = "FRR")
    mc <- mc %>% select(GEOID, COUNTY_MEAN, REGION_MEAN, DIFF, Z)
    bs <- ifelse(mc$Z > 2, "Bright spot", NA)
    bs <- ifelse(mc$Z < -2, "Dark spot", bs)
    bs[is.na(bs)] <- "Average"
    mc$BS <- bs
    colnames(mc) <- c("GEOID", 
                      paste0(toupper(crops[i]), "_CM"),
                      paste0(toupper(crops[i]), "_RM"), 
                      paste0(toupper(crops[i]), "_DIFF"),
                      paste0(toupper(crops[i]), "_Z"),
                      "geometry",
                      paste0(toupper(crops[i]), "_BS"))
    
    mc <- as(mc, "Spatial")

    rgdal::writeOGR(mc, paste0("./out/shp/", crops[i], ".shp"), driver = "ESRI Shapefile", layer = "diff")
    
    cty <- merge(cty, mc, by = "GEOID", all = T)
    
  }
  
  return(cty)
  
  if (save == T) {
    cty <- as(cty, "Spatial")
    rgdal::writeOGR(cty, "./out/diff.shp", driver = "ESRI Shapefile", layer = "Diff")
  }
}

compute_nat_diff <- function(model, data, level) {
  
  r <- model
  b0 <- r$summary.fixed[1,1]
  
  # extract area-specific residuals
  n <- length(unique(data[,"ID"]))
  asr <- r$summary.random[["ID"]][1:n, c(1, 2)]  # E_ij = u_0jk + v_0jk
  colnames(asr) <- c("ID", "COUNTY_MEAN")
  
  m <- merge(county_sub, asr, by = "ID", all=T)
  m$NATIONAL_MEAN <- mean(data$YIELD, na.rm=T)

  m$DIFF <- m$COUNTY_MEAN - m$NATIONAL_MEAN
  m$Z <- NA
  
  lrr_ids <- unique(data$FRR)
  for (i in 1:length(lrr_ids)) {
    id <- lrr_ids[[i]]
    df <- m %>% filter(FRR == id)
    df$Z <- (df$DIFF - mean(df$DIFF, na.rm=T))/sd(df$DIFF, na.rm=T)
    m$Z[m$FRR == id] <- df$Z
  }
  
  return(m)
  
}

# functinos for attribute analysis

prep_attr_data <- function(crop, scale = T) {
  
  # prepare data
  train <- prep_data(crop)[[1]]
  county_sub <- prep_data(crop)[[2]]
  null <-  readRDS(paste0("./out/", crop, ".RDS"))
  modr <- readRDS(paste0("./out/models/mf_", crop, ".RDS"))
  
  # crop sqkm
  null$PERC_CROP <- null[,paste0(toupper(crop), "_SQKM")]/null$TOTAL_SQKM * 100
  null$PERC_AG <- null$AG_SQKM/null$TOTAL_SQKM * 100
  perc <- null %>% select(GEOID, YEAR, PERC_CROP, PERC_AG)
  
  # load Census data
  hf <- readRDS("./data/BSDS_standatt_operated_v3.RDS")
  hf$YEAR <- hf$year
  hfm <- merge(hf, perc, by = c("GEOID", "YEAR"), all = T) %>%
    filter(YEAR %in% c(2007, 2012, 2017)) %>% # only for years in train dataset
    filter(GEOID %in% unique(train$GEOID)) # and GEOIDs in dataset
  hfm <- hfm %>%
    select(-c(year, YEAR, county)) %>% # remove null model attr
    group_by(GEOID) %>% 
    summarize_all(funs(mean), na.rm=T)
  diff <- compute_diff(modr, train, "FRR")
  regdf <- merge(hfm, diff, by = "GEOID", all = T)
  # df <- regdf %>% select(-c(cons_wet_acres, manure_acres, male, insur_op, labor_n, 
  #                          acres_per_op, exp, insur_acres,
  #                            cty_cl, cty_clp, cty_al, cty_pe, ph_corn, perc_cl, perc_clp,
  #                            perc_p, perc_pe, perc_o, perc_wle, perc_t)) 
  df <- regdf %>% select(-c(cty_cl, cty_clp, cty_al, cty_pe, ph_corn, perc_cl, perc_clp,
                            perc_p, perc_pe, perc_o, perc_wle, perc_t)) 

  frag <- readRDS("C:/Users/eburchf/OneDrive - Emory University/WF/RR_LSI/LSI/out/indices/LSM_AREA_MN_AG.RDS") %>%
    group_by(GEOID) %>% summarize(AREA_MN_AG = mean(LSM_AREA_MN_AG, na.rm=T))
  sdi <- readRDS("C:/Users/eburchf/OneDrive - Emory University/WF/RR_LSI/LSI/out/indices/LSM_SHDI_AG.RDS") %>%
    group_by(GEOID) %>% summarize(SDI = mean(LSM_SHDI_AG, na.rm=T))
  ed <- readRDS("C:/Users/eburchf/OneDrive - Emory University/WF/RR_LSI/LSI/out/indices/LSM_ED_AG.RDS") %>%
    group_by(GEOID) %>% summarize(ED = mean(LSM_ED_AG, na.rm=T))
  lpi <- readRDS("C:/Users/eburchf/OneDrive - Emory University/WF/RR_LSI/LSI/out/indices/LSM_LPI_AG.RDS") %>%
    group_by(GEOID) %>% summarize(LPI_AG = mean(LSM_LPI_AG, na.rm=T))
  pnc <- readRDS("C:/Users/eburchf/OneDrive - Emory University/WF/RR_LSI/LSI/out/indices/PNC.RDS") %>%
    group_by(GEOID) %>% summarize(PNC = mean(PNC, na.rm=T))
  df <- merge(df, frag, by = "GEOID", all=T)
  df <- merge(df, sdi, by = "GEOID", all = T)
  df <- merge(df, lpi, by = "GEOID", all = T)
  df <- merge(df, pnc, by = "GEOID", all = T)
  df <- merge(df, ed, by = "GEOID", all = T)
  st_geometry(df) <- df$geometry
  
  df <- df %>% select(-c("ID", "REGION_MEAN", "COUNTY_MEAN", "STATEFP"))
  
  # df <- df %>% select(-c(PERC_AG, comm_sales, income,
  #                        crop_sales, herb_acres, insect_acres, fert_acres, part_owner,
  #                        irrig,
  #                        income_farm,
  #                        cons_wetlands,
  #                        tenant,
  #                        chem # too collinear
  # )) %>% na.omit()
  # 
  # scale data for comparison across models 
  
  if (scale == T) {
    scaled <- df %>% select(-c(GEOID, DIFF, Z, FRR, geometry)) 
    st_geometry(scaled) <- NULL
    scaled <- as.data.frame(scale(scaled)) 
    scaled <- cbind.data.frame(df %>% select(c(GEOID, DIFF, Z, FRR, geometry)), scaled)
    out <- scaled
  } else {
    out <- df
  }
  return(out)
}

plot_coef_lm <- function(model) {
  
  coef <- tidy(summary(model))
  coef <- cbind.data.frame(coef$term, coef$estimate, coef$std.error)
  colnames(coef) <- c("VARIABLE", "MEAN", "SE")
  coef <- coef %>% filter(VARIABLE != "(Intercept)")
  
  ggplot(coef) +
    geom_errorbarh(height = 0, aes(xmin = MEAN - SE, 
                                   xmax = MEAN + SE, 
                                   y = reorder(VARIABLE, desc(MEAN)))) +
    project_theme +
    geom_point(aes(x = MEAN, y = reorder(VARIABLE, desc(MEAN))), size = 2) +
    geom_vline(xintercept = 0) +
    labs(x = "",
         y = "")
}

run_linear_model <- function(scaled){
  
  ## 75% of the sample size
  smp_size <- floor(0.75 * nrow(scaled))
  
  ## set the seed to make your partition reproducible
  set.seed(123)
  train_ind <- sample(seq_len(nrow(scaled)), size = smp_size)
  
  train <- scaled[train_ind, ]
  test <- scaled[-train_ind, ]
  
  out <- lm(DIFF ~ ., 
            data = train %>% select(-c(geometry, GEOID, FRR, Z)))
  preds <- out %>% predict(test)
  diag <- data.frame(RMSE = RMSE(preds, test$DIFF),
                     R2 = caret::R2(preds, test$DIFF))
  
  out <- lm(DIFF ~ ., 
            data = scaled %>% select(-c(geometry, GEOID, FRR, Z)))
  
  out <- list(out, diag)
  return(out)
  
  
}
