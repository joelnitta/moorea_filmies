# Set up captions

# Follow Journal of Plant Research (JPR) style:
# - All figures are to be numbered using Arabic numerals.
# - Figures should always be cited in text in consecutive numerical order.
# - Figure parts should be denoted by lowercase letters (a, b, c, etc.).
# - If an appendix appears in your article and it contains one or more figures, 
#   continue the consecutive numbering of the main text. 
# - Do not number the appendix figures, "A1, A2, A3, etc." 
# - Figures in online appendices (Electronic Supplementary Material) should, 
#   however, be numbered separately

# - Figures
figure_full <- captioner::captioner(prefix = "Fig.", suffix = "")
figure_full(name = "filmy-habit", caption = "filmy-habit")
figure_full(name = "filmy-dt", caption = "filmy-dt")

# - Tables
table_full <- captioner::captioner(prefix = "Table")

# - SI figures
s_figure_full <- captioner::captioner(prefix = "Fig. S", suffix = "")

# - SI tables
s_table_full <- captioner::captioner(prefix = "S", auto_space = FALSE, suffix = " Table. ")

# Make short versions of citation functions (just the number)
figure <- pryr::partial(figure_full, display = "cite")
table <- pryr::partial(table_full, display = "cite")
s_figure <- pryr::partial(s_figure_full, display = "cite")
s_table <- pryr::partial(s_table_full, display = "cite")
