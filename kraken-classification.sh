#!/bin/bash

# Set up variables
DB_NAME=dash-cam-db
DB_DIR=./kraken2-data

rm -rf $DB_DIR

# Create the database directory
mkdir -p $DB_DIR/$DB_NAME

kraken2/kraken2-inst/kraken2-build --db $DB_DIR/$DB_NAME --download-taxonomy

# Build the custom Kraken2 database
kraken2/kraken2-inst/kraken2-build --add-to-library data/influenza.fna --db $DB_DIR/$DB_NAME
kraken2/kraken2-inst/kraken2-build --add-to-library data/lassa.fna --db $DB_DIR/$DB_NAME
kraken2/kraken2-inst/kraken2-build --add-to-library data/measles.fna --db $DB_DIR/$DB_NAME
kraken2/kraken2-inst/kraken2-build --add-to-library data/sars-cov-2.fna --db $DB_DIR/$DB_NAME
kraken2/kraken2-inst/kraken2-build --add-to-library data/rotavirus.fna --db $DB_DIR/$DB_NAME

# Build the Kraken2 database index
kraken2/kraken2-inst/kraken2-build --build --db $DB_DIR/$DB_NAME
