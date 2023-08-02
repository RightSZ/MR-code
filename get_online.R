while(TRUE){
  message_to_next <<- TRUE
  error_to_next <<- FALSE
  try({withCallingHandlers(exp <- expr, 
                           message = function(c) if (stringr::str_detect(as.character(c),"Failed to")) message_to_next <<- FALSE)
    error_to_next <<- TRUE})
  if(message_to_next == TRUE&error_to_next == TRUE) { break }
}

