#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyalert)
library(shinyjs)
library(shinyBS)
library(shinycssloaders)
library(ggplot2)
library(dplyr)
library(readr)
library(bslib)


thesis <- read_csv("thesis.csv")
pubs <- read_csv("pubs.csv")
teach <- read_csv("teach.csv")
projects <- read_csv("projects.csv")
projects <- projects %>% mutate(
  logo = c('<img src="https://pbs.twimg.com/profile_images/1675084495626682374/Kl6AOIau_400x400.png" height="30"></img>',
           '<img src="https://mlt.org/wp-content/uploads/2020/11/Biogen-Logo.png" height="50"></img>',
           '<img src="https://bigtencrc.org/wp-content/uploads/2021/05/uicc.png" height="30"></img>',
           '<img src="https://assets-global.website-files.com/63bc6cdba9784b6ec05f51aa/6473c746453365df38678d84_63176094cf736012a6a4f0ee_surest%25E2%2584%25A2_logo_lockup_purp_RGB.png" height="30"></img>',
           '<img src="https://upload.wikimedia.org/wikipedia/commons/thumb/c/c3/Mayo_Clinic_logo.svg/1200px-Mayo_Clinic_logo.svg.png" height="50"></img>'
           )
  )


ui <- shinyUI(dashboardPage(
  skin = "midnight",
  #controlbar = dashboardControlbar(collapsed = FALSE, skinSelector()),
  # header and title page
  dashboardHeader(title = "Zhaoliang Zhou's Personal Website", titleWidth = 400, disable = FALSE,
                  tags$li(class = "dropdown", tags$a(href = "https://www.linkedin.com/in/zhaoliang-zhou-646b52115/", 
                                                     icon("linkedin"), 
                                                     "LinkedIn", 
                                                     target = "_blank")),
                  tags$li(class = "dropdown", tags$a(href = "https://github.com/zhaoliang-zhou", 
                                                     icon("github"), 
                                                     "GitHub", 
                                                     target = "_blank"))
                  #target options 
                  #_self - Default. Opens the document in the same window/tab as it was clicked
                  #_blank - Opens the document in a new window or tab
                  #_parent - Opens the document in the parent frame
                  #_top - Opens the document in the full body of the window
  ),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem(text = "About", tabName = "about", icon = icon("bullseye")), 
      menuItem(text = "Research", tabName = "research", icon = icon("cube"),  
               menuSubItem("Papers", tabName = "pub"),
               menuSubItem("Projects", tabName = "proj")),
      menuItem(text = "Teaching", tabName = "teach", icon = icon("users")),
      menuItem(text = "Misc.", tabName = "misc", icon = icon("heart"),
               menuSubItem("Favorite movies/TV", tabName = "favmov")),
      menuItem(text = "CV", tabName = "cv", icon = icon("code"))
    )
  ),
  
  dashboardBody(
    tabItems(
      tabItem(tabName = "about", uiOutput("indexpage")),
      
      tabItem(tabName = "pub", uiOutput("pubpage")),
      
      tabItem(tabName = "proj", DT::DTOutput("projtable")),
      
      tabItem(tabName = "teach", 
              DT::DTOutput("teachtable")
              ),
      
      # misc UI 
      tabItem(tabName = "favmov", 
              fluidRow(
                infoBoxOutput("BB", width = 3),
                infoBoxOutput("moneyheist", width = 3),
                infoBoxOutput("GOT", width = 3),
                infoBoxOutput("band", width = 3)),
              fluidRow( # create a new row
                infoBoxOutput("shaw", width = 3),
                infoBoxOutput("last", width = 3))
        ),
      
      
      # CV UI
      tabItem(tabName = "cv", tags$iframe(style = "height: 1100px; width: 100%; scrolling = yes",
                                          src = "CV-ZZ.pdf"))
    )
    
  )
  
))


server <- function(input, output, session) {
  
  # About page
  output$indexpage <- renderUI({
    sidebarLayout(position = "left",
                  mainPanel(width = 12,
                               tags$img(src = "umngrad.jpg", height="30%", width="30%", 
                                        style="display: block; margin-left: auto; margin-right: auto;"),
                               h4(tags$b("Current Ph.D student in Biostatistics, University of Illinois - Chicago (UIC)"),
                                  align="center",
                                  style = "color:white"), 
                               
                               hr(),
                               
                               h5(tags$b("Education"),style = "color:white"),
                               h5("Ph.D in Biostatistics, University of Illinois - Chicago (UIC), 2023 - present"),
                               h5("M.S. in Statistics, University of Minnesota - Twin Cities, 2020 - 2022"),
                               h5("B.A. with Economics major and Statistics and Data Science minor, St. Olaf College, 2016 - 2020"),
                               
                               hr(),
                               
                               h5(tags$b("Experiences"),style = "color:white"),
                               h5("Commodity Analyst, AbbVie, 09/2023 - 05/2024"),
                               h5("Biostatistics Intern, Biogen, 06/2023 - 08/2023"),
                               h5("Research Assistant, UIC, 08/2022 - 05/2023"),
                               h5("Health Equity Analytic Intern, Bind Benefits (Surest Health Plan), 06/2021 - 08/2021"),
                               h5("Biostatistics Intern, Mayo, 05/2019 - 08/2019"),
                               
                               hr(),
                               
                               h5(tags$b("Contact"),style = "color:white"),
                               h5("Email: zz81@uic.edu")
                  #end of main panel             
                  ),
                 mainPanel(
                  # must have a mainPanel() here, otherwise error message: "mainPanel" missing
                  # h4("Hi, welcome to my website!")
                  )
                
    )
  # end of About page  
  })
  
  # publication page
  output$pubpage <- renderUI({
    tabsetPanel(
      tabPanel(tags$b("Thesis"), DT::DTOutput("thesistable")),
      tabPanel(tags$b("Publications"), DT::DTOutput("pubstable"))
    )
  })
  
  output$thesistable <- DT::renderDT({
    DT::datatable(thesis,style = 'bootstrap4')
  })
  
  output$pubstable <- DT::renderDT({
    DT::datatable(pubs,style = 'bootstrap4', escape = FALSE)
  })

  output$projtable <- DT::renderDT({
    DT::datatable(projects,style = 'bootstrap4', escape = FALSE)
  })
  
  output$teachtable <- DT::renderDT({
    DT::datatable(teach,style = 'bootstrap4', escape = FALSE)
  })
  
  # misc favorite TV series
  output$BB <- renderInfoBox({
    infoBox("Breaking Bad", "", 
            div(img(src = "breakingbad.jpg", width = 200), 
                style = "text-align: center;"),
            icon = icon("thumbs-up"),
            href = "https://www.youtube.com/watch?v=VaOt6tXyf2Y")
  })
  output$moneyheist <- renderInfoBox({
    infoBox("Money Heist", "", 
            div(img(src = "moneyheist.jpg", width = 200), 
                style = "text-align: center;"),
            icon = icon("thumbs-up"),
            href = "https://www.youtube.com/watch?v=_InqQJRqGW4")
  })
  output$GOT <- renderInfoBox({
    infoBox("Game of Thrones", "", 
            div(img(src = "got.jpg", width = 150), 
                style = "text-align: center;"),
            icon = icon("thumbs-up"),
            href = "https://www.youtube.com/watch?v=KPLWWIOCOOQ")
  })
  output$band <- renderInfoBox({
    infoBox("Band of Brothers", "", 
            div(img(src = "band.jpg", width = 200), 
                style = "text-align: center;"),
            icon = icon("thumbs-up"),
            href = "https://www.youtube.com/watch?v=KKRBAFlN5ww")
  })
  output$shaw <- renderInfoBox({
    infoBox("Shawshank's Redemption", "", 
            div(img(src = "shaw.jpg", width = 200), 
                style = "text-align: center;"),
            icon = icon("thumbs-up"),
            href = "https://www.youtube.com/watch?v=PLl99DlL6b4")
  })
  output$last <- renderInfoBox({
    infoBox("The last of us", "", 
            div(img(src = "last.jpg", width = 200), 
                style = "text-align: center;"),
            icon = icon("thumbs-up"),
            href = "https://www.youtube.com/watch?v=uLtkt8BonwM")
  })
# end of server  
}

# Run the application 
shinyApp(ui = ui, server = server)




# notes
# list of icons avaliable: https://fontawesome.com/icons
# for files uploaded from local, make sure create a new file called www and put the files into that folder
# this website used some of the ideas/structures from https://github.com/Xiaozhu-Zhang1998/website-using-shiny/tree/main 
# other package for themes: {fresh} https://unleash-shiny.rinterface.com/beautify-with-fresh#beautify-with-fresh

