#ffmpeg -f image2 -framerate 30 -i lj-%03d.png -vcodec libx264 -crf 22 -pix_fmt yuv420p video_4K.mp4
ffmpeg -f image2 -framerate 30 -i pov/lj-%03d.png -vcodec libx264 -crf 22 -pix_fmt yuv420p lj_all_inside_HD.mp4

