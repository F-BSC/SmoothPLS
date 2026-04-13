This is a resubmission. In this version I have addressed the comments from Uwe Ligges:
* Removed the redundant phrase starting the Description field and started with "Performs...".
* Added single quotes around acronyms and software names (e.g., 'SmoothPLS') or expanded them to resolve spellcheck notes.
* Added the methodological reference (Aguilera et al., 2010) with its DOI directly in the Description field.

## Test environments
* local Windows 11 install, R 4.4.1
* win-builder (devel)

## R CMD check results
0 errors | 0 warnings | 1 note

* NOTE: checking CRAN incoming feasibility
  * This is a new submission.

## URL checks
* NOTE: `urlchecker::url_check()` reports a 403 Forbidden error for 
`https://www.decathlon.fr/`. This is a false positive. The website uses 
strict anti-bot protections (like Cloudflare) that block automated HTTP HEAD 
requests from R. 
This URL is perfectly valid and accessible via a standard web browser.
