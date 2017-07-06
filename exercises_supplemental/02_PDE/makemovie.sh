mencoder \
	"mf://*.png" -mf fps=15 \
	-o ausgabevideo.avi \
	-ovc lavc \
	-lavcopts vcodec=mpeg4:vbitrate=8000
