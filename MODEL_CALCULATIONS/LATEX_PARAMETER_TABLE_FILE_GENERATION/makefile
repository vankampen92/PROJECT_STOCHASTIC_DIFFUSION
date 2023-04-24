FILE_NAME = Parameter_Table_Master

latex: 
	latex ${FILE_NAME} 	

bib:
	bibtex ${FILE_NAME}

ps:
	dvips -o ${FILE_NAME}.ps ${FILE_NAME}
	gv ${FILE_NAME}.ps

pdf: 
	pdflatex ${FILE_NAME}
	acroread ${FILE_NAME}.pdf

ps2pdf: 
	ps2pdf ${FILE_NAME}.ps
	acroread ${FILE_NAME}.pdf

clean:
	rm *~ *dvi *log *bbl *blg *aux
