commit:
	git commit -a --allow-empty-message -m ''
	git push

doc:
	doxygen

clean:
	rm -rf Documenation
