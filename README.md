# Sepsis-3 criteria and epidemiology in AmsterdamUMCdb

This repository contains files for implementing the Sepsis-3 definition in the freely-accessible Amsterdam University Medical Centers Database. The folder `sepsis3_gigabyte' contains code accompanying [our manuscript](http://dx.doi.org/10.46471/gigabyte.45) published in the journal Gigabyte. The folder `sepsis3_epidemiology' contains complementary code for a forthcoming manuscript describing sepsis epidemiology in AmsterdamUMCdb.

You must first have access to AmsterdamUMCdb, including the data. Further instructions for this can be found on the [Amsterdam Medical Data Science](https://amsterdammedicaldatascience.nl/) website or [AmsterdamUMCdb GitHub](https://github.com/AmsterdamUMC/AmsterdamUMCdb). To run code in this repositor, refer to either of the folders `sepsis3_gigabyte' or `sepsis3_epidemiology'.

The database is primarily in Dutch. It is generally easier to work with itemids rather than the (Dutch) item names, in any table. The database contains a lookup dictionary of all terms, instructions for installing the `amsterdamumcdb` package that contains this is given [here](https://github.com/AmsterdamUMC/AmsterdamUMCdb/tree/master/setup-amsterdamumcdb). It may be worth familiarising yourself with some common terms (e.g. Neurochirurgie = Neurosurgery, Vaatchirurgie = Vascular Surgery) ahead of time.

We encourage other researchers (or anyone interested in this code from a data or clinical perspective) to submit an [issue](https://github.com/tedinburgh/sepsis3-amsterdamumcdb/issues) on this repository if they would like to report or fix any bugs or to contribute improvements.
