library(shiny)

source("r.value.r")

shinyServer(function(input, output){
  
  makeTable <- reactive(
  {
    if (is.null(input$file$name)) stop("Choose .csv file")
    if (is.na(input$m)) stop("enter m : number of hypotheses in primary study")
    if (is.na(input$c2)) stop("value for c2 is missing")
    if (is.na(input$l00)) stop("value for l_00 is missing")
    if((1<=input$q)|(input$q<=0)) stop("q should be in the interval (0,1)")
    if (input$pChk)
      if (input$var=="use.t")
        if (is.na(input$t)) stop("value for t is missing")
    if (!is.na(input$t))
    {
      if ((1<=input$t)|(input$t<=0))
        stop("t should be in the interval (0,1)")
      else
        tt <- input$t
    }
    else
      tt <- NULL
    
    input_table <- read.csv(file = input$file$datapath, header=T, stringsAsFactors =F,
                            colClasses = c("character","numeric","numeric"))
    
    params <- list(m = input$m, k = nrow(input_table), l00 = input$l00, c2 = input$c2,
                   q = input$q, variation = input$var, t = tt)
    
    if (input$pChk)
      rv <- r.value(p1 = input_table[,2], p2 = input_table[,3], m = input$m, c2 = input$c2,
                    l00 = input$l00, variation = input$var, tt = tt , Q = input$q)
    else
      rv <- r.value(p1 = input_table[,2], p2 = input_table[,3], m = input$m, c2 = 0.5,
                    l00 = input$l00)
    
    output_table <- cbind(input_table,rv,ifelse(rv < input$q,"*",""))
    colnames(output_table) <- c("Index\\Name","Primary p-value","Secondary p-value",
                                "r-value",paste("Significant (r-value <=",ifelse(input$pChk,input$q,0.05),")"))
    return(list(table=output_table,params=params))
  })
  
  output$log <- renderText({
    if (input$goButton == 0)
      return()
    else
    {
      params <- isolate({makeTable()$params})
      if (isolate(input$pChk))
      {
        out <- paste('Data input: <em>m</em> = ',params[[1]],'; <em>R1</em> = ',params[[2]],'<br>
                 Parameters: <em>l_00</em> = ',params[[3]],'; <em>c2</em> = ',params[[4]],'; <em>q</em> = ',params[[5]],'; variation: ',params[[6]],ifelse(params[[6]]=="use.t",paste(" ; <em>t</em> = ",params[[7]],sep=''),""),sep='')
        
        if (params[[6]] != "use.t" & !is.null(params[[7]])) 
          out <- paste0(out,"<div class='warn'>Warning: Threshold t is ignored.</div>")
        if (params[[6]] == "use.t" & !is.null(params[[7]]))
        {
          if (params[[7]] <= (1-params[[4]])/(1-params[[3]]*(1-params[[4]]*params[[5]])) * params[[5]]/params[[1]] )
            out <- paste0(out,"<div class='warn'>Warning: since t < c(q)q/m, no modification to the original r-value computation was necessary (see section Derivation and Properties in the article)</div>")
          if (params[[7]] >= (1-params[[4]])/(1-params[[3]]*(1-params[[4]]*params[[5]])) * params[[5]]/(1+sum(1/(1:(params[[1]]-1)))))
            out <- paste0(out,"<div class='warn'>Warning: for the selected threshold t, the 'use.t' variation won't lead\nto more discoveries than the 'use.m.star' variation.</div>")
        } 
      }
      else out <- paste('Data input: <em>m</em> = ',params[[1]],'; <em>R1</em> = ',params[[2]],'<br>
                    Parameters: <em>l_00</em> = ',params[[3]],'; <em>c2</em> = ',params[[4]],sep='')
      return(out)
    }
  })
  
output$table <- renderTable({
  if (input$goButton == 0)
    return()
  else
  {
    table <- isolate({makeTable()$table})
    return(table)
  }
},display = c("s","s","G","G","G","s"),align=c("r","r","r","r","r","c"),include.rownames=F
 ,sanitize.colnames.function=function(x) gsub("<=","&#8804",x))

output$downloadtable <- downloadHandler(
  filename = function()
  {
    if (input$goButton == 0)
      return()
    else
    {
      params <- isolate({makeTable()$params})
      if (isolate(input$pChk)) paste('RepOut_m.',params[1],'_R1.',params[2],'_l00.',params[3],'_c2.',params[4],'_q.',params[5],'_var.',params[6],ifelse(params[6]=="use.t",paste("_t.",params[7],sep=''),""),'_time.', Sys.time(),'.csv', sep='')
      else                     paste('RepOut_m.',params[1],'_R1.',params[2],'_l00.',params[3],'_c2.',params[4],"_time.",Sys.time(),'.csv', sep='')
    }
  }
    ,
  content = function(file) {
    if (input$goButton == 0)
      return()
    else
    {
      table <- isolate({makeTable()$table})
      write.csv(table, file, row.names=F)
    }
  })
  outputOptions(output, "downloadtable", suspendWhenHidden=FALSE)
})