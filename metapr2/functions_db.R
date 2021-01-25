# Database libraries
  library(DBI)
  library(RMySQL)

# =================================================
  
db_info <- function(database_name,
                    file_cnf = "mysql/my.cnf")  {
# To do, try https://theautomatic.net/2019/06/25/how-to-hide-a-password-in-r-with-the-keyring-package/


  db_pr2_google <- list(dbname='pr2',
                        default.file=file_cnf,
                        groups="google-pr2")
  
  db_metapr2_google <- list(dbname='metapr2',
                        default.file=file_cnf,
                        groups="google-metapr2")
  dbinfo <- NULL

  if (database_name == "pr2_google") {db_info=db_pr2_google}
  if (database_name == "metapr2_google") {db_info=db_metapr2_google}

    return(db_info)
}

# ====================================

db_connect <- function(db_info)  {

db <- dbConnect(MySQL(),  default.file=db_info$default.file,
                          groups=db_info$groups,
                          dbname=db_info$dbname)
return(db)
}


# ====================================

db_disconnect <- function(db_con)  {
    dbDisconnect(db_con)
    return(TRUE)
}
