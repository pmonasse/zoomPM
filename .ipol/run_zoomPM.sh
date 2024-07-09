#!/bin/bash
# Arguments: zoom color
# -zoom: zoom factor
# -color: handle color image?

zoom=$1
color=$2
color=${color/true/-c}
color=${color/false/}

echo "$bin/build/zoomPM -z $zoom $color input_0.png zoom.png" > demo_failure.txt
#$bin/build/zoomPM -z $zoom $color input_0.png zoom.png
