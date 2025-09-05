# this script creates an app demonstrating page_navbar

library(shiny)
library(bslib)
library(tidyverse)
library(AnophelesModel)
library(cowplot)

#setwd(dirname(rstudioapi::getSourceEditorContext()$path))
dirInputs="Inputs"

source("./LLIN_cascade_inputs.R")

source("./cascade_functions.R")
source("./compute_vectorialCapacity.R")
source("./wrapper_cascade_functions.R")

all_results<- read.csv(file.path(dirInputs,"fitted_parameters_posteriormax.csv"))
eht_description<- read.csv(file.path(dirInputs,"EHT_description.csv"))
names(eht_description)=c("Experimental hut trial", "Location","Year","Mosquito species", "New generation ITNs", "Durability" )
activity_patterns_Ref<- readRDS(file.path(dirInputs,"activity_patterns_withRefs_small.rds"))


activity_noRhythms =def_activity_patterns()
activity_noRhythms$humans_in_bed=rep(1, length(activity_noRhythms$HBI))
activity_noRhythms$humans_indoors=rep(1, length(activity_noRhythms$HBI))


host_table = as.data.frame(def_host_params())

ui <- page_navbar(

  
  
  
  title = span(img(src = "logosmall.png"), "ITN Effectiveness cascades"),
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
          title = "Description",
          p("This dashboard provides effectiveness cascades for two insecticide-treated nets, namely Interceptor G2 and Olyset Plus. 
            Effectiveness is measured in terms of vectorial capacity reduction. 
            The cascade enables the user to quantify the loss of effectiveness attributable to five factors, namely entomological efficacy, usage at distribution, functional survival, insecticidal durability and in-bed exposure (as defined in Definition tab).
            The user can use the panel on the left-hand side to adapt the assumptions and tailor the results to their setting of interest."),
          
          p(strong("Interpreting the outcome:"),"Using the default assumptions, the effectiveness of Interceptor G2 nets was above 90% when only entomological efficacy was considered. It dropped below 60% when all other dimensions were included. In this example, the main factors responsible for the decline in effectiveness were in-bed exposure (responsible for a drop of 13 points), insecticidal durability (responsible for a drop of 12 points) and imperfect usage at distribution (responsible for a drop of 10 points)"),
          p("A short summary of the assumptions and definitions is provided in the Assumption tab and all methodological details are available in the following manuscript"),
          p( strong("Cascades of effectiveness of new-generation insecticide-treated nets against malaria, from entomological trials to real-life conditions")),
          p("Clara Champagne, Jeanne Lemant, Alphonce Assenga, Ummi A. Kibondo, Ruth G. Lekundayo, Emmanuel Mbuba, Jason Moore, Joseph B. Muganga, Watson S. Ntabaliba, Olukayode G. Odufuwa, Johnson Kyeba Swai, Maria Alexa, Roland Goers, Monica Golumbeanu, Nakul Chitnis, Amanda Ross, Raphael N'Guessan, Sarah Moore, Emilie Pothin"),
          p( strong("Contact:") , "clara.champagne@swisstph.ch")
        ), 
        
        nav_panel(
           title = "ITN cascades",
          layout_columns(
             card(card_header("Interceptor G2"),plotOutput(outputId = "IG2")),
             card(card_header("Olyset Plus"),plotOutput(outputId = "OlysetPlus"))
        )
         ),
        nav_panel(
          title = "Assumptions and definitions",
          navset_card_tab(
            nav_panel(
              title = "Assumptions",
              p("These dashboards are associated with the following manuscript ", strong("'Quantifying the effectiveness of new generation bed nets against malaria, from entomological trials to real life conditions'"), " to which the user should refer for methodological details.
            The sources and assumptions underlying the default parameter values are as follows:",
                p(strong("Entomological efficacy:"), "Capacity of a given vector control product to reduce mosquito transmission potential, through the alteration of mosquito biting, mortality and/or fecundity, under ideal operational conditions. It can be quantified as the relative reduction in vectorial capacity through deploying the vector control product under ideal operational conditions."),
                p(strong("ITN effectiveness:"), "Capacity of a given ITN product to reduce mosquito transmission potential under user conditions. It differs from entomological efficacy in that insecticidal durability, functional survival, ITN use and in-bed exposure are also accounted for. It can be quantified as the relative reduction in vectorial capacity through deploying the ITN under realistic operational conditions."),
                p(strong("Vectorial capacity:"), "Total number of potentially infectious bites originating from all the mosquitoes biting a single perfectly infectious (i.e. all mosquito bites result in infection) human on a single day. It represents the potential for a given mosquito population to transmit malaria."),
                
                p(strong("EHT:"), "Experimental hut trial dataset used to estimate entomological efficacy. Characteristics of the different EHT included are displayed on the subsequent tab."),
                p(strong("ITN use:"),"Proportion of the population sleeping under a net the previous night. Default values are estimated with usage data from" ,
                  tags$a("Mosha et al. 2022", href="https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(21)02499-5/fulltext", target="_blank" ), " and ",
                  tags$a("Mosha et al. 2024", href="https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(23)00420-6/fulltext", target="_blank" )),
                p(strong("Functional survival:"), "An estimate of the physical lifespan of an ITN product in the field. It is measured as an ITN that is present in use and in serviceable condition.
                It is the result of the ITN's physical integrity declining due to damage and attrition, that results in nets being discarded when the user no longer regards them as useful.
                  It is represented with a Weibull function with user-specific half-life and shape parameters. Default values are estimated with usage data from" ,
                  tags$a("Mosha et al. 2022", href="https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(21)02499-5/fulltext", target="_blank" ), " and ",
                  tags$a("Mosha et al. 2024", href="https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(23)00420-6/fulltext" , target="_blank"),
                  ".
                  In practice, for other settings, online tools (e.g.", tags$a("ITN Quantification", href="https://allianceformalariaprevention.com/itn-quantification/", target="_blank" ) ,") are available to estimate functional survival and ITN needs over time using survey data and stock-and-flow models (",
                  tags$a("Bertozzi-Villa et al. 2021", href="https://www.nature.com/articles/s41467-021-23707-7", target="_blank" ), 
                  tags$a("Koenker et al. 2023", href="https://malariajournal.biomedcentral.com/articles/10.1186/s12936-023-04609-z", target="_blank" ),")."),
                p(strong("Anopheles bionomics:"), "Parous rate, human blood index, sac rate and endophagy in the model correspond to the values for the chosen species (default is Anopheles gambiae). 
              Values are extracted from ",
                  tags$a("Golumbeanu et al. 2024", href="https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011609", target="_blank"),
                  "via",
                  tags$a("AnophelesModel", href="https://github.com/SwissTPH/AnophelesModel", target="_blank"), "and",
                  tags$a("AnophelesBionomics", href="https://github.com/SwissTPH/AnophelesBionomics", target="_blank")),
                p(strong("Insecticidal durability:"), "Retention of insecticidal efficacy of the ITN. It is estimated as the difference in entomological efficacy between unwashed/new ITNs and 20-times washed/field-aged ITNs in EHTs."),
                p(strong("In-bed exposure:"), "'The proportion of vector bites occurring indoors during sleeping hours, for an unprotected individual, which represents the maximum possible personal protection any intervention targeting sleeping spaces could provide' (",
                  tags$a("Monroe et al. 2020", href="https://doi.org/10.1186/s12936-020-03271-z", target="_blank" ),"). It is calculated by combining activity patterns for humans and mosquitoes as in",
                  tags$a("Golumbeanu et al. 2024", href="https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011609", target="_blank")),
                p(strong("Activity patterns:"), "Using activity rhythm data averaged for the chosen country (default: Tanzania). For vector data, all available species and surveys are averaged. For human data, all available surveys are averaged.
              Values are extracted from",
                  tags$a("Golumbeanu et al. 2024", href="https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011609", target="_blank"),
                  "via",
                  tags$a("AnophelesModel", href="https://github.com/SwissTPH/AnophelesModel", target="_blank"), ". The data is displayed on the subsequent tab."),
            )
            ),
            nav_panel(
              title = "Experimental hut trial data",
              #plotOutput(outputId = "plot_rhythms")
              dataTableOutput(outputId = "EHTtable")
            ),
            nav_panel(
            title = "Activity patterns",
            plotOutput(outputId = "plot_rhythms")
            )
        
        )
        
        )
    )
   )
  
  
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
    updateNumericInput(session, "attrition_shape_ig2", value = 2.4)
    updateSelectInput(session, "EHT_ig2", selected = "Martin et al., An. gambiae")
  })
  
  observeEvent(input$reset_default_pbo, {
    updateSelectInput(session, "species_pbo", selected = "Anopheles gambiae")
    updateSelectInput(session, "rhythms_hu_pbo", selected = "Tanzania")
    updateSelectInput(session, "rhythms_mosq_pbo", selected = "Tanzania")
    updateNumericInput(session, "usage_pbo", value = 0.75)
    updateNumericInput(session, "attrition_pbo", value = 1.7)
    updateNumericInput(session, "attrition_shape_pbo", value = 1.9)
    updateSelectInput(session, "EHT_pbo", selected = "Martin et al., An. gambiae")
  })
  
  observeEvent(input$set_pbo_as_IG2, {
    updateSelectInput(session, "species_pbo", selected = input$species_ig2)
    updateSelectInput(session, "rhythms_hu_pbo", selected = input$rhythms_hu_ig2)
    updateSelectInput(session, "rhythms_mosq_pbo", selected = input$rhythms_mosq_ig2)
    updateNumericInput(session, "usage_pbo", value = input$usage_ig2)
    updateNumericInput(session, "attrition_pbo", value = input$attrition_ig2)
    updateNumericInput(session, "attrition_shape_pbo", value = input$attrition_shape_ig2)
  })
  
  observeEvent(input$set_IG2_as_pbo, {
    updateSelectInput(session, "species_ig2", selected = input$species_pbo)
    updateSelectInput(session, "rhythms_hu_ig2", selected = input$rhythms_hu_pbo)
    updateSelectInput(session, "rhythms_mosq_ig2", selected = input$rhythms_mosq_pbo)
    updateNumericInput(session, "usage_ig2", value = input$usage_pbo)
    updateNumericInput(session, "attrition_ig2", value = input$attrition_pbo)
    updateNumericInput(session, "attrition_shape_ig2", value = input$attrition_shape_pbo)
  })
  
    is_pbo_notas_ig2 <- eventReactive(eventExpr = input$update_pbo,
                                     ignoreNULL = F,
                                     {
                                       return(input$rhythms_hu_ig2 != input$rhythms_hu_pbo | input$rhythms_mosq_ig2!=input$rhythms_mosq_pbo )                         
                                       
                                     })
  
  is_ig2_notas_pbo <- eventReactive(eventExpr = input$update_ig2,
                                 ignoreNULL = F,
                                 {
                                   return(input$rhythms_hu_ig2 != input$rhythms_hu_pbo | input$rhythms_mosq_ig2!=input$rhythms_mosq_pbo )                         
                                   
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
                             activity_ig2=activity_function(select_rythms_ig2()$my_rhythms_mean,species=input$species_ig2)                          
                             
                             outputs_posteriormax = all_results[(all_results$reference==input$EHT_ig2 & all_results$netType=="Interceptor G2"),]
                             VCreduc_all=run_anophelesModelforCascade(ent_params=ent_params_ig2,
                                                                      host_params=host_params_ig2,
                                                          outputs_posteriormax=outputs_posteriormax,
                                                          activity_noRhythms=activity_noRhythms,usage=input$usage_ig2, exposure=activity_ig2,
                                                          attrition=input$attrition_ig2, kappa_attrition=input$attrition_shape_ig2)        
                             
                             #print(VCreduc_all)
                             return(VCreduc_all)
                           })
  
  
  output$IG2 <- renderPlot({
      
    plot_bars_effectiveness_LLIN(df_VCred=cascade_ig2(),colorfinal="orange")#+
      # theme(axis.text = element_text(size = 12) )
     
  #  }

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
                                activity_pbo=activity_function(select_rythms_pbo()$my_rhythms_mean,species=input$species_pbo)                          
                                
                                outputs_posteriormax = all_results[(all_results$reference==input$EHT_pbo & all_results$netType=="Olyset Plus"),]
                                
                                VCreduc_all=run_anophelesModelforCascade(ent_params=ent_params_pbo,
                                                                         host_params=host_params_pbo,
                                                                         outputs_posteriormax=outputs_posteriormax,
                                                                         activity_noRhythms=activity_noRhythms,usage=input$usage_pbo, exposure=activity_pbo,
                                                                         attrition=input$attrition_pbo, kappa_attrition=input$attrition_shape_pbo)        
                                
                                return(VCreduc_all)
                              })
 
  output$OlysetPlus <- renderPlot({
    

    plot_bars_effectiveness_LLIN(df_VCred=cascade_pbo(), 
                                 colorfinal="dodgerblue")#+
      #theme(axis.text = element_text(size = 12) )
    
    #  }
    
  })

  
  ################
  # Activity patterns
  
  
  output$plot_rhythms <- renderPlot({
    
    if(is_pbo_notas_ig2() | is_ig2_notas_pbo()){
      plot_grid(
        rhythms_function_plot(select_rythms_ig2()$my_rhythms, select_rythms_ig2()$my_rhythms_mean)+
          guides(color=guide_legend(nrow=3)),
        rhythms_function_plot(select_rythms_pbo()$my_rhythms, select_rythms_pbo()$my_rhythms_mean)+
          guides(color=guide_legend(nrow=3)),
        labels=c("Interceptor G2","Olyset Plus"), scale = 0.9
      )
    } else{
      rhythms_function_plot(select_rythms_ig2()$my_rhythms, select_rythms_ig2()$my_rhythms_mean)
    }
    
  })

  output$EHTtable <-renderDataTable(eht_description,
                                    options = list(searching = FALSE))
  
}

shinyApp(ui, server)