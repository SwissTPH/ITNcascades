# this script creates an app demonstrating page_navbar

library(shiny)
library(bslib)
library(tidyverse)
library(AnophelesModel)

#setwd(dirname(rstudioapi::getSourceEditorContext()$path))
dirInputs="Inputs"

source("./LLIN_cascade_inputs.R")

source("./cascade_functions.R")
source("./compute_vectorialCapacity.R")
source("./wrapper_cascade_functions.R")

results_kibondo_w<- readRDS(file.path(dirInputs,"output_kibondo_washed.rds"))
results_kibondo_unw<- readRDS(file.path(dirInputs,"output_kibondo_unwashed.rds"))

results_bit055_w<- readRDS(file.path(dirInputs,"output_bit055_washed.rds"))
results_bit055_unw<- readRDS(file.path(dirInputs,"output_bit055_unwashed.rds"))

results_bit103_unw<- readRDS(file.path(dirInputs,"output_bit103_unwashed.rds"))
results_bit103_w<- readRDS(file.path(dirInputs,"output_bit103_washed.rds"))

results_bit059_w<- readRDS(file.path(dirInputs,"output_bit059_washed.rds"))
results_bit059_unw<- readRDS(file.path(dirInputs,"output_bit059_unwashed.rds"))

results_bit080_unw<- readRDS(file.path(dirInputs,"output_bit080_unwashed.rds"))
results_bit080_w<- readRDS(file.path(dirInputs,"output_bit080_washed.rds"))

activity_patterns_Ref<- readRDS(file.path(dirInputs,"activity_patterns_withRefs_small.rds"))




activity_noRhythms =def_activity_patterns()
activity_noRhythms$humans_in_bed=rep(1, length(activity_noRhythms$HBI))
activity_noRhythms$humans_indoors=rep(1, length(activity_noRhythms$HBI))


host_table = as.data.frame(def_host_params())

ui <- page_navbar(

  
  
  
  title = span(img(src = "logosmall.png"), "LLIN Effectiveness cascades"),
  #title = "LLIN Effectiveness cascades",
  # nav_panel(
  #   title = "Introduction", 
  #   #tags$h2("Heading"),
  #   tags$p("Modelling estimates of various interventions against", tags$i("Plasmodium vivax")  ,"in Lao PDR.")
  #   ),
  # 
  #############################
  # LLIN cascades
  #############################

    layout_sidebar(
      sidebar = sideBar_InputsCascade,
      navset_card_tab(
         nav_panel(
           title = "LLIN cascades",
          layout_columns(
             card(card_header("Interceptor G2"),plotOutput(outputId = "IG2")),
             card(card_header("OlysetPlus"),plotOutput(outputId = "OlysetPlus"))
        )
         ),
        nav_panel(
          title = "Assumptions",
          navset_card_tab(
            nav_panel(
              title = "Assumptions",
              p("These dashboards are associated with the following manuscript ", strong("'Quantifying the effectiveness of new generation bed nets against malaria, from entomological trials to real life conditions'"), " to which the user should refer for methodological details.
            The sources and assumptions underlying the parameter values are as follows:",
                p(strong("Bednet use:"),"Proportion of the population sleeping under a bednet. Default values are estimated with usage data from" ,
                  tags$a("Mosha et al. 2022", href="https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(21)02499-5/fulltext" ), " and ",
                  tags$a("Mosha et al. 2024", href="https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(23)00420-6/fulltext" )),
                p(strong("Anopheles species:"), "Using bionomics (parous rate, human blood index, sac rate, endophily and endophagy) corresponding to the given species. 
              Values are extracted from",
                  tags$a("Golumbeanu et al. 2024", href="https://www.biorxiv.org/content/10.1101/2023.10.17.562838v1"),
                  "via",
                  tags$a("AnophelesModel", href="https://github.com/SwissTPH/AnophelesModel")),
                p(strong("Activity patterns:"), "Using activity rhythm data averaged for the chosen country (default: Tanzania). For vector data, all available species and surveys are averaged. For human data, all available surveys are averaged.
              Values are extracted from",
                  tags$a("Golumbeanu et al. 2024", href="https://www.biorxiv.org/content/10.1101/2023.10.17.562838v1"),
                  "via",
                  tags$a("AnophelesModel", href="https://github.com/SwissTPH/AnophelesModel"), "The data is displayed on the subsequent tabs."),
                p(strong("EHT:"), "Experimental hut trial dataset used to estimate entomological efficacy"),
                p(strong("Attrition  (half life, in years):"), "Time until which less than half of the originally covered population are still using a bednet (assuming a Weibull function). Default values are estimated with usage data from" ,
                  tags$a("Mosha et al. 2022", href="https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(21)02499-5/fulltext" ), " and ",
                  tags$a("Mosha et al. 2024", href="https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(23)00420-6/fulltext" )),
                p(strong("Attrition  (shape):"), "Shape parameter of the Weibull distribution specified above. Default values are estimated with usage data from" ,
                  tags$a("Mosha et al. 2022", href="https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(21)02499-5/fulltext" ), " and ",
                  tags$a("Mosha et al. 2024", href="https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(23)00420-6/fulltext" )),
            )
            ),
            nav_panel(
              title = "Activity patterns IG2",
              plotOutput(outputId = "plot_rhythms_ig2")
            ),
            nav_panel(
              title = "Activity patterns OlysetPlus",
              plotOutput(outputId = "plot_rhythms_pbo")
            )
          )#

        ),
        nav_spacer(),
        nav_panel(
            title = "About",
            p("This dashboard is related to the following manuscript"),
            p( strong("Cascades of effectiveness of new-generation insecticide-treated nets against malaria, from entomological trials to real-life conditions")),
p("Clara Champagne, Jeanne Lemant, Alphonce Assenga, Ummi A. Kibondo, Ruth G. Lekundayo, Emmanuel Mbuba, Jason Moore, Joseph B. Muganga, Watson S. Ntabaliba, Olukayode G. Odufuwa, Johnson Kyeba Swai, Maria Alexa, Roland Goers, Monica Golumbeanu, Nakul Chitnis, Amanda Ross, Sarah Moore, Emilie Pothin"),
p( strong("Contact:") , "clara.champagne@swisstph.ch")
          )
        )
        
    
    )
  #  )
  
  
)

server <- function(input, output, session) {
  
  ################
  # Cascade reset
  observeEvent(input$reset_default_ig2, {
    updateSelectInput(session, "species_ig2", selected = "Anopheles gambiae")
    updateSelectInput(session, "rhythms_hu_ig2", selected = "Tanzania")
    updateSelectInput(session, "rhythms_mosq_ig2", selected = "Tanzania")
    updateNumericInput(session, "usage_ig2", value = 0.69)
    updateNumericInput(session, "attrition_ig2", value = 2.4)
    updateNumericInput(session, "attrition_shape_ig2", value = 1.9)
    updateSelectInput(session, "EHT_ig2", selected = "BIT103")
    

  })
  
  observeEvent(input$reset_default_pbo, {
    updateSelectInput(session, "species_pbo", selected = "Anopheles gambiae")
    updateSelectInput(session, "rhythms_hu_pbo", selected = "Tanzania")
    updateSelectInput(session, "rhythms_mosq_pbo", selected = "Tanzania")
    updateNumericInput(session, "usage_pbo", value = 0.75)
    updateNumericInput(session, "attrition_pbo", value = 1.7)
    updateNumericInput(session, "attrition_shape_pbo", value = 0.8)
    updateSelectInput(session, "EHT_pbo", selected = "BIT103")
    
    
    
  })
  ################
  # Cascade IG2
 
  
 select_rythms_ig2 <- eventReactive(eventExpr = input$update_ig2,
                                    ignoreNULL = F,
                                    {
                                      return(rhythms_function(country_hu=input$rhythms_hu_ig2, country_mosq = input$rhythms_mosq_ig2,
                                                              my_activity_patterns=activity_patterns_Ref))                         
                                      
                                    })
  
cascade_ig2 <- eventReactive(eventExpr = input$update_ig2,
                           ignoreNULL = F,
                           {
                             
                             ent_params_ig2 = def_vector_params(mosquito_species = input$species_ig2, verbose = FALSE)
                             
                             host_table_ig2=host_table
                             host_table_ig2$species_name = input$species_ig2
                             host_params_ig2 = def_host_params(mosquito_species = input$species_ig2, vec_params = ent_params_ig2, host_table = host_table_ig2, verbose = FALSE)
                             #activity_ig2=rhythms_function(country_hu=input$rhythms_hu_ig2, country_mosq = input$rhythms_mosq_ig2)                          
                             activity_ig2=activity_function(select_rythms_ig2()$my_rhythms_mean)                          
                             
                             
                             if(input$EHT_ig2=="Kibondo"){
                               my_results_unw=results_kibondo_unw
                               my_results_w=results_kibondo_w
                               my_insecticide_id=2
                             } else if (input$EHT_ig2=="BIT103"){
                               my_results_unw=results_bit103_unw
                               my_results_w=results_bit103_w
                               my_insecticide_id=3
                             } else if (input$EHT_ig2=="BIT080"){
                               my_results_unw=results_bit080_unw
                               my_results_w=results_bit080_w
                               my_insecticide_id=2
                             }
                             
                             VCreduc_all=run_anophelesModelforCascade(ent_params=ent_params_ig2,
                                                                      host_params=host_params_ig2,
                                                          activity_p=activity_ig2 , activity_noRhythms=activity_noRhythms,
                                                          results=my_results_unw, results_washed=my_results_w, 
                                                          usage=input$usage_ig2, 
                                                          insecticide_id=my_insecticide_id, 
                                                          attrition=input$attrition_ig2, kappa_attrition=input$attrition_shape_ig2)        
                             
                             return(VCreduc_all)
                           })
  
  
  output$IG2 <- renderPlot({
      
    plot_bars_effectiveness_LLIN(df_VCred=cascade_ig2(),colorfinal="orange")#+
      # theme(axis.text = element_text(size = 12) )
     
  #  }

  })
  
  
  output$plot_rhythms_ig2 <- renderPlot({
    
    rhythms_function_plot(select_rythms_ig2()$my_rhythms, select_rythms_ig2()$my_rhythms_mean)
    
  })
  ################
  # Cascade OlysetPlus
 
  select_rythms_pbo <- eventReactive(eventExpr = input$update_pbo,
                                     ignoreNULL = F,
                                     {
                                       return(rhythms_function(country_hu=input$rhythms_hu_pbo, country_mosq = input$rhythms_mosq_pbo,
                                                               my_activity_patterns=activity_patterns_Ref) )                         
                                       
                                     })
 cascade_pbo <- eventReactive(eventExpr = input$update_pbo,
                              ignoreNULL = F,
                              {
                                
                                ent_params_pbo = def_vector_params(mosquito_species = input$species_pbo)
                                host_table_pbo=host_table
                                host_table_pbo$species_name = input$species_pbo
                                host_params_pbo = def_host_params(mosquito_species = input$species_pbo, vec_params = ent_params_pbo, host_table = host_table_pbo, verbose = FALSE)
                                #activity_pbo=rhythms_function(country_hu=input$rhythms_hu_pbo, country_mosq = input$rhythms_mosq_pbo)                          
                                activity_pbo=activity_function(select_rythms_pbo()$my_rhythms_mean)                          
                                
                                
                                if(input$EHT_pbo=="BIT055"){
                                  my_results_unw=results_bit055_unw
                                  my_results_w=results_bit055_w
                                  my_insecticide_id=2
                                } else if(input$EHT_pbo=="Odufuwa"){
                                  my_results_unw=results_bit059_unw
                                  my_results_w=results_bit059_w
                                  my_insecticide_id=2
                                } else if(input$EHT_pbo=="BIT103"){
                                  my_results_unw=results_bit103_unw
                                  my_results_w=results_bit103_w
                                  my_insecticide_id=2
                                }
                                
                                VCreduc_all=run_anophelesModelforCascade(ent_params=ent_params_pbo,
                                                                         host_params=host_params_pbo,
                                                                         activity_p=activity_pbo , activity_noRhythms=activity_noRhythms,
                                                                         results=my_results_unw, results_washed=my_results_w, 
                                                                         usage=input$usage_pbo, 
                                                                         insecticide_id=my_insecticide_id, 
                                                                         attrition=input$attrition_pbo, kappa_attrition=input$attrition_shape_pbo)        
                                
                                return(VCreduc_all)
                              })
 
  output$OlysetPlus <- renderPlot({
    

    plot_bars_effectiveness_LLIN(df_VCred=cascade_pbo(), 
                                 colorfinal="dodgerblue")#+
      #theme(axis.text = element_text(size = 12) )
    
    #  }
    
  })
  
  output$plot_rhythms_pbo <- renderPlot({
    
    rhythms_function_plot(select_rythms_pbo()$my_rhythms, select_rythms_pbo()$my_rhythms_mean)
    
  })

  
}

shinyApp(ui, server)