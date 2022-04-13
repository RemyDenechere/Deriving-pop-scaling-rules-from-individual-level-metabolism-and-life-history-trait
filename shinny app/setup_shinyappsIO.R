#install.packages('rsconnect')

#rsconnect::setAccountInfo(name='remydenechere-centreforoceanlife',
#                          token='B75C1478A3D75E80DB4604D07A263A82',
#                          secret='+tCURePW6xe7AmxSpIMh4adPwM4OyVqylDmLY58e')
library(rsconnect)
# upload rmax app
rsconnect::deployApp('C:/Users/rden/OneDrive - Danmarks Tekniske Universitet/Ph.D/Project_1/Deriving-pop-scaling-rules-from-individual-level-metabolism-and-life-history-trait/shinny app/Rmax')
