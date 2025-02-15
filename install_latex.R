# Install latex packages using tinytex

# This would happen automatically anyways when rendering the Rmd to pdf, 
# but requires downloading packages and may not work during updates of Tex Live. 
# Better to install to the docker image once and keep them there.

tinytex::tlmgr_update()

latex_packages <- c(
  "amsmath",
  "auxhook",
  "bigintcalc",
  "bitset",
  "booktabs",
  "caption",
  "colortbl",
  "environ",
  "etexcmds",
  "etoolbox",
  "fancyvrb",
  "filehook",
  "float",
  "framed",
  "geometry",
  "gettitlestring",
  "grffile",
  "hycolor",
  "hyperref",
  "infwarerr",
  "intcalc",
  "kvdefinekeys",
  "kvoptions",
  "kvsetkeys",
  "latex-amsmath-dev",
  "letltxmacro",
  "lm-math",
  "ltxcmds",
  "makecell",
  "mdwtools",
  "multirow",
  "pdfescape",
  "pdflscape",
  "pdftexcmds",
  "refcount",
  "rerunfilecheck",
  "setspace",
  "siunitx",
  "stringenc",
  "tabu",
  "threeparttable",
  "threeparttablex",
  "trimspaces",
  "ucharcat",
  "ulem",
  "unicode-math",
  "uniquecounter",
  "varwidth",
  "wrapfig",
  "xcolor",
  "zapfding")

tinytex::tlmgr_install(latex_packages)