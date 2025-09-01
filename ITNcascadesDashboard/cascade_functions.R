plot_bars_effectiveness_LLIN      <- function(df_VCred, colorfinal="dodgerblue"){
  
  
  # Create a dataframe with the right format
  
  bars_df=data.frame(variable="bites", value=100) %>% 
    rbind(df_VCred) %>%
    rbind(df_VCred %>% filter(variable=="EfficacyAttritionUsageDecayRhythms")%>% mutate(variable="effectiveness"))%>%
    mutate(
      background_value = lag(value),
      #percent          = 100 * (background_value-value) / background_value
      percent          = 100 * (value) / background_value
    )
  
  bars_df$value[bars_df$variable %in% c("effectiveness")]=NA
  bars_df$variable=factor(bars_df$variable, levels=c("bites", "Efficacy","EfficacyUsage",
                                                     "EfficacyAttritionUsage","EfficacyAttritionUsageDecay","EfficacyAttritionUsageDecayRhythms",
                                                     "effectiveness"),
                          labels=c("No\nvectors",
                                   "Entomological\nefficacy",
                                   "Usage\nat distribution",
                                   "Functional\nsurvival",
                                   "Insecticidal\ndurability",
                                   "In-bed\nexposure",
                                   "Effectiveness"))
  bars_df$point_difference=round(bars_df$background_value-bars_df$value)
  
  bars_background <- ggplot(bars_df) +
    geom_col(mapping = aes(variable ,
                           background_value,
                           fill = (variable %in% c("Effectiveness"))
    ),
    color = "black", position = "dodge"
    ) +
    scale_fill_manual(values = c("grey47", colorfinal)) +
    guides(fill = "none") +
    geom_text(data = subset(bars_df, variable %in% c("Effectiveness")),
              mapping = aes(
                variable, background_value,
                label = background_value %>% round()%>% paste0("%")
              ),
              position = position_stack(vjust = 0.15),
    ) +
    scale_y_continuous(label=function(x) paste0(x, "%")) +
    labs(x="",
         y="Reduction in vectorial capacity") +
    theme_classic()
  
  my_plot=bars_background +
    geom_col(mapping = aes(variable, value),
             fill = "lightgrey",
             color= "black") +
    geom_text(data = bars_df %>% filter(variable !="No\nvectors" & variable != "Effectiveness"),
              mapping = aes(
                variable, value,
                #label = percent %>% round() %>% paste0("%\nless")
                label = value %>% round() %>% paste0("%")
              ),
              position = position_stack(vjust = 0.15),
    ) +
    geom_label(data = bars_df %>% filter(variable !="No\nvector" & variable != "Effectiveness"),
               mapping = aes(
                 variable,
                 label = paste0("-",point_difference,"pts"),
                 y=value+5
               ), #nudge_x = 0.5,#color="white",
               position = position_stack()
    )
  return(my_plot)
}
