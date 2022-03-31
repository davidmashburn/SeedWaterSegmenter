mkdir SeedWaterSegmenter.iconset
sips -z 16 16     SeedWaterSegmenter1024Circle.png --out SeedWaterSegmenter.iconset/icon_16x16.png
sips -z 32 32     SeedWaterSegmenter1024Circle.png --out SeedWaterSegmenter.iconset/icon_16x16@2x.png
sips -z 32 32     SeedWaterSegmenter1024Circle.png --out SeedWaterSegmenter.iconset/icon_32x32.png
sips -z 64 64     SeedWaterSegmenter1024Circle.png --out SeedWaterSegmenter.iconset/icon_32x32@2x.png
sips -z 128 128   SeedWaterSegmenter1024Circle.png --out SeedWaterSegmenter.iconset/icon_128x128.png
sips -z 256 256   SeedWaterSegmenter1024Circle.png --out SeedWaterSegmenter.iconset/icon_128x128@2x.png
sips -z 256 256   SeedWaterSegmenter1024Circle.png --out SeedWaterSegmenter.iconset/icon_256x256.png
sips -z 512 512   SeedWaterSegmenter1024Circle.png --out SeedWaterSegmenter.iconset/icon_256x256@2x.png
sips -z 512 512   SeedWaterSegmenter1024Circle.png --out SeedWaterSegmenter.iconset/icon_512x512.png
cp SeedWaterSegmenter1024Circle.png SeedWaterSegmenter.iconset/icon_512x512@2x.png
iconutil -c icns SeedWaterSegmenter.iconset
rm -R SeedWaterSegmenter.iconset
