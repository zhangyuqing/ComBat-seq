rm(list=ls())
setwd("~/yuqingz/ComBat_seq")
sapply(c("ballgown"), require, character.only=TRUE)

load("geuvadisbg.rda")
geuvadisbg@expr
#Error: no slot of name "expr" for this object of class "ballgown"
