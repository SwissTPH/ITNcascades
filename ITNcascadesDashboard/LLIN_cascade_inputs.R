sideBar_InputsCascade=sidebar(
  title="Parameters",
  p("For definitions and assumptions on parameter values please see Assumptions tab"),
  accordion(
    open = FALSE,
    accordion_panel(
      title = "Interceptor G2",
      actionButton( inputId = "update_ig2",label = "Update plot"
      ),
      numericInput( inputId="usage_ig2",  label="ITN use",  value=0.69,  min = 0,  max = 1,
                    step = 0.05,
                    width = NULL
      ),
      selectInput( inputId="species_ig2",  label="Anopheles bionomics",
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
                  choices = c("Kibondo et al.", "Assenga et al., Tanzania","Assenga et al., Côte d'Ivoire", "BIT080", 
                              "Martin et al., An. gambiae", "Martin et al., An. funestus",
                              "Sovegnon et al.", "Nguessan et al." ), selected = "Martin et al., An. gambiae"),
      
      numericInput( inputId="attrition_ig2",  label="Functional survival (half life, in years)",  value=2.4,  min = 0,  max = 5,
                    step = 0.1,
                    width = NULL
      ),
      numericInput( inputId="attrition_shape_ig2",  label="Functional survival (shape)",  value=2.4,  min = 0,  max = 5,
                    step = 0.1,
                    width = NULL
      ),
      actionButton("reset_default_ig2", "Default values IG2"),
      actionButton("set_IG2_as_pbo", "Same values as Olyset Plus")
      
    ),
    accordion_panel(
      title = "OlysetPlus",
      actionButton( inputId = "update_pbo",label = "Update plot"
      ),
      numericInput( inputId="usage_pbo",  label="ITN use",  value=0.75,  min = 0,  max = 1,
                    step = 0.05,
                    width = NULL
      ),
      selectInput( inputId="species_pbo",  label="Anopheles bionomics",
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
                  choices = c( "Assenga et al., Tanzania","Assenga et al., Côte d'Ivoire", "BIT055", "Odufuwa et al.",
                              "Martin et al., An. gambiae", "Martin et al., An. funestus"), selected = "Martin et al., An. gambiae"),
      
      numericInput( inputId="attrition_pbo",  label="Functional survival (half life, in years)",  value=1.7,  min = 0,  max = 5,
                    step = 0.1,
                    width = NULL
      ),
      numericInput( inputId="attrition_shape_pbo",  label="Functional survival (shape)",  value=1.9,  min = 0,  max = 5,
                    step = 0.1,
                    width = NULL
      ),
      actionButton("reset_default_pbo", "Default values Olyset Plus"),
      actionButton("set_pbo_as_IG2", "Same values as IG2")
    ),
    
  ),
  
  p(""),

)
