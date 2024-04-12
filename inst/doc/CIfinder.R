## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----global, include=FALSE----------------------------------------------------
library(CIfinder)
library(kableExtra)

## ----eval = FALSE-------------------------------------------------------------
#  install.packages("CIfinder")

## ----echo=FALSE---------------------------------------------------------------
bdat <- data.frame(Case = c(31, 3, 34), Control = c(12, 32, 44))
rownames(bdat) <- c("Poor_Signature", "Good_Signature", "Total")
kable(bdat, caption = "Table 1. Breast cancer data from a molecular signature classifier") %>%
  kable_classic(full_width = F, html_font = "Cambria")

## ----gn-----------------------------------------------------------------------
library(CIfinder)
ppv_npv_ci(x1 = 31, n1 = 34, x0 = 32, n0 = 44, prevalence = 0.07,
           method = "gart and nam")

## ----zhou1--------------------------------------------------------------------
ppv_npv_ci(x1 = 31, n1 = 34, x0 = 32, n0 = 44, prevalence = 0.07,
           method = "zhou")

## ----zhou2--------------------------------------------------------------------
ppv_npv_ci(x1 = 31, n1 = 34, x0 = 32, n0 = 44, prevalence = 0.07,
           method = "zhou",
           continuity.correction = TRUE)

## ----boot---------------------------------------------------------------------
ppv_npv_ci(x1 = 31, n1 = 34, x0 = 32, n0 = 44, prevalence = 0.07,
           method = "boot")

## ----fieller------------------------------------------------------------------
ppv_npv_ci(x1 = 31, n1 = 34, x0 = 32, n0 = 44, prevalence = 0.07,
           method = "fieller")

## ----echo=FALSE---------------------------------------------------------------
bdat <- data.frame(Case = c(31, 3, 34), Control = c(0, 44, 44))
rownames(bdat) <- c("Poor_Signature", "Good_Signature", "Total")
kable(bdat, caption = "Table 1. Example of a special testing data") %>%
  kable_classic(full_width = F, html_font = "Cambria")

## ----gn2----------------------------------------------------------------------
ppv_npv_ci(x1 = 31, n1 = 34, x0 = 44, n0 = 44, prevalence = 0.07,
           method = "gart and nam")

## ----walter-------------------------------------------------------------------
ppv_npv_ci(x1 = 31, n1 = 34, x0 = 44, n0 = 44, prevalence = 0.07,
           method = "walter")

## ----zhou3--------------------------------------------------------------------
ppv_npv_ci(x1 = 31, n1 = 34, x0 = 44, n0 = 44, prevalence = 0.07,
           method = "zhou",
           continuity.correction = TRUE)

