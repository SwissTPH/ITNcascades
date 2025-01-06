sideBar_InputsCascade=sidebar(
  title="Parameters",
  accordion(
    open = FALSE,
    accordion_panel(
      title = "Interceptor G2",
      actionButton( inputId = "update_ig2",label = "Update plot"
      ),
      numericInput( inputId="usage_ig2",  label="Bednet use",  value=0.69,  min = 0,  max = 1,
                    step = 0.05,
                    width = NULL
      ),
      selectInput( inputId="species_ig2",  label="Anopheles species",
                   choices = AnophelesModel::list_all_species(), 
                   #choices = c("Anopheles gambiae", "Anopheles funestus", "Anopheles arabiensis"),#AnophelesModel::list_all_species(), 
                   selected = "Anopheles gambiae"),

      selectInput( inputId="rhythms_hu_ig2",  label="Human activity patterns",  
                   choices = c("Cote D'Ivoire", "Haiti", "Kenya",
                               "Mozambique",  "Tanzania","Zambia"), selected = "Tanzania"),
      
      selectInput( inputId="rhythms_mosq_ig2",  label="Mosquito activity patterns",  
                  choices = c("Benin", "Burkina Faso","Cameroon", "Democratic Republic of Congo", "Equatorial Guinea", 
                              "Erithrea", "Gabon", "Ghana", "Haiti", "India", "Kenya",
                              "Mozambique", "Nigeria", "PNG", "Sao Tome", "Senegal",
                              "Tanzania", "Thailand", "The Gambia", "Uganda","Zambia", "Zimbabwe"), selected = "Tanzania"),
      
      selectInput("EHT_ig2", label = "EHT",
                  choices = c("Kibondo", "BIT103", "BIT080"), selected = "BIT103"),
      
      numericInput( inputId="attrition_ig2",  label="Attrition  (half life, in years)",  value=2.4,  min = 0,  max = 5,
                    step = 0.1,
                    width = NULL
      ),
      numericInput( inputId="attrition_shape_ig2",  label="Attrition  (shape)",  value=2.4,  min = 0,  max = 5,
                    step = 0.1,
                    width = NULL
      ),
      actionButton("reset_default_ig2", "Default values IG2")
      
    ),
    accordion_panel(
      title = "OlysetPlus",
      actionButton( inputId = "update_pbo",label = "Update plot"
      ),
      numericInput( inputId="usage_pbo",  label="Bednet use",  value=0.75,  min = 0,  max = 1,
                    step = 0.05,
                    width = NULL
      ),
      selectInput( inputId="species_pbo",  label="Anopheles species",
                   choices = AnophelesModel::list_all_species(), 
                   selected = "Anopheles gambiae"),
      
      selectInput( inputId="rhythms_hu_pbo",  label="Human activity patterns",  
                   choices = c("Cote D'Ivoire", "Haiti", "Kenya",
                               "Mozambique",  "Tanzania","Zambia"), selected = "Tanzania"),
      
      selectInput( inputId="rhythms_mosq_pbo",  label="Mosquito activity patterns",  
                   choices = c("Benin", "Burkina Faso","Cameroon", "Democratic Republic of Congo", "Equatorial Guinea", 
                               "Erithrea", "Gabon", "Ghana", "Haiti", "India", "Kenya",
                               "Mozambique", "Nigeria", "PNG", "Sao Tome", "Senegal",
                               "Tanzania", "Thailand", "The Gambia", "Uganda","Zambia", "Zimbabwe"), selected = "Tanzania"),
      selectInput("EHT_pbo", label = "EHT",
                  choices = c("BIT055", "BIT103", "Odufuwa"), selected = "BIT103"),
      
      numericInput( inputId="attrition_pbo",  label="Attrition (half life, in years)",  value=1.7,  min = 0,  max = 5,
                    step = 0.1,
                    width = NULL
      ),
      numericInput( inputId="attrition_shape_pbo",  label="Attrition (shape)",  value=1.9,  min = 0,  max = 5,
                    step = 0.1,
                    width = NULL
      ),
      actionButton("reset_default_pbo", "Default values PBO")
    ),
    
  ),
  
  p(""),

)
