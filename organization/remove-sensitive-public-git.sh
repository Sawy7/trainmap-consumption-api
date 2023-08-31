#!/bin/sh
# Remove sensitive data from git
git-filter-repo --invert-paths --path code-doc/olo-opava.geojson
git-filter-repo --invert-paths --path testing-data
