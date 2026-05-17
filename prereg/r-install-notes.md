# R Package Installation Notes

## Summary
Successfully installed 97 out of the required R packages. Missing system dependencies prevented installation of some packages.

## Failed Packages & Required System Dependencies

### 1. curl (7.1.0)
- **Error**: `libcurl` not found
- **System dependency**: `libcurl4-openssl-dev` (Debian/Ubuntu)
- **Command**: `sudo apt-get install -y libcurl4-openssl-dev`
- **Affected dependents**: httr, rvest, kableExtra, lintr, nhanesA, xml2, timeROC

### 2. xml2 (1.5.2)
- **Error**: `libxml-2.0` not found
- **System dependency**: `libxml2-dev` (Debian/Ubuntu)
- **Command**: `sudo apt-get install -y libxml2-dev`
- **Affected dependents**: rvest, kableExtra, lintr, nhanesA, timeROC

### 3. rms (*) - Wildcard version
- **Error**: Package marked with asterisk version; failed to download
- **Note**: rms is typically available on CRAN but the version resolution failed
- **Workaround**: May require specifying explicit version or checking package availability
- **Affected dependents**: riskRegression, pec, timeROC

### 4. riskRegression (2026.03.11)
- **Error**: Dependency 'rms' not available
- **Note**: Cannot install until rms is resolved
- **Requires**: rms package to be available

### 5. timeROC (0.4.1)
- **Error**: Dependency 'pec' failed (which requires rms and riskRegression)
- **Note**: Indirect dependency on rms through pec package

### 6. kableExtra (1.4.0)
- **Error**: Dependency 'xml2' installation failed
- **Note**: Only needed for table formatting in reports; non-critical

### 7. lintr (3.3.0-1)
- **Error**: Dependency 'xml2' installation failed
- **Note**: Linting tool; non-critical for analysis

### 8. nhanesA (1.4.1)
- **Error**: Dependencies 'rvest' and 'xml2' failed
- **Note**: Critical for data acquisition from NHANES
- **Priority**: High - install system deps to get this working

## Successfully Installed (97 packages)

survey, cmprsk, survival, pROC, dcurves, mice, rpart, rpart.plot, dplyr, tidyr, purrr, readr, ggplot2, targets, tarchetypes, testthat, and 81 dependencies.

## Installation Path Forward

1. **Immediate** (requires sudo password): Install system dependencies
   ```bash
   sudo apt-get install -y libcurl4-openssl-dev libxml2-dev libgsl-dev
   ```

2. **After system deps**: Re-run renv::install() for failed packages
   ```r
   renv::install(c("nhanesA", "riskRegression", "timeROC", "kableExtra", "lintr", "arrow"))
   ```

3. **Note on arrow**: Large package (~400MB pre-compiled); consider installing separately if time is limited

## Critical vs Optional

- **Critical for analysis**: nhanesA (data loading), riskRegression (competing risks), timeROC (time-dependent ROC)
- **Quality tools** (optional): kableExtra, lintr
- **Performance** (optional): arrow (for large data I/O)

## Environment Notes

- R version: 4.3.3
- renv version: 1.2.3
- Ubuntu/Debian system libraries needed for compilation from source
- All available packages are cached; re-installation will be fast once system deps are installed
