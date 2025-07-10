# CDS bootstrapping library:
Library implementing 3 methods to bootstrap survival probability curve and intensity curve from market CDS spreads from different maturities:
- Model with no accrual ;
- Model with accrual ;
- Jarrow Turnbull approximation.

The library works the following:
LIB_BOOTSTRAP_CDS.m is like a demo which calls bootstrapCDS.m and given the used flag (1, 2 or 3), it calls either bootstrapCDS_NOaccrual.m, bootstrapCDS_accrual.m or bootstrap_JT.m.
 ## NB: to discount the future expected cashflow, the discount curve is needed. Therefore in LIB_BOOTSTRAP_CDS.m we start by bootstrapping the discount curve using bootstrap.m, readExcelData.m and the market data given in an excel file (or .mat).
