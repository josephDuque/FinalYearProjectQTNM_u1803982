default : git 

.PHONY: git
git:
	git add .
	git commit -m "$m"
	git push -u origin main

.PHONY: ammendGit
ammendGit:
	git reset --soft "HEAD^"
	git commit --amend
	git push -u origin main

