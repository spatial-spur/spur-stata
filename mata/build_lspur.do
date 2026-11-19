* Run in a fresh Stata session from the spur-stata repository root.
* Optional argument: directory for the compiled library.
version 14
args output_dir
if `"`output_dir'"' == "" local output_dir "."
set more off
mata: mata clear
mata: mata set matastrict on

local sources : dir "mata" files "*.mata"
foreach source of local sources {
    quietly do "mata/`source'"
}
mata: mata mlib create lspur, dir(`"`output_dir'"') size(64) replace
mata: mata mlib add lspur *(), dir(`"`output_dir'"')
mata: mata mlib index
display "SPUR MATA LIBRARY BUILT"
