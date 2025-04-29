library(httr)
library(jsonlite)
library(xml2)
 
server <- "http://localhost:8000}"
ext <- "/api/get_sequence/PP760209"
 
r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
 
stop_for_status(r)
 
# use this if you get a simple nested list back, otherwise inspect its structure
# head(data.frame(t(sapply(content(r),c))))
head(fromJSON(toJSON(content(r))))