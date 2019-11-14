import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import matplotlib.patches as patches
import sys

file = '/home/mkirsche/phase_block_lengths_enc004.txt'
ofn = 'phaseblocks_enc004.png'

if len(sys.argv) == 3:
  file = sys.argv[1]
  ofn = sys.argv[2]

with open(file) as f:
  lines = f.readlines()
  chrCount = len(lines)/2
  print('Number of chromosomes: ' + str(chrCount))

  pixelsPerInch = 80
  totalHeightInches = 8
  totalWidthInches = 12
  totalPaddingPixels = 30
  textWidthPixels = 100
  totalHeight = pixelsPerInch * totalHeightInches
  totalWidth = pixelsPerInch * totalWidthInches + textWidthPixels
  propHeightOccupied = 0.8
  totalOccupied = int(totalHeight * propHeightOccupied)
  chrHeight = int(totalOccupied / chrCount)
  spaceHeight = int((totalHeight - totalOccupied) / (chrCount + 1))

  # Iterate over the chromosomes to find the longest one for scaling purposes
  maxChrLength = 0
  for i in range(0, len(lines), 2):
    next = lines[i+1].split()
    nextInt = [int(val) for val in next]
    totalLength = max(nextInt)
    if totalLength > maxChrLength:
      maxChrLength = totalLength
  print('Maximum chromosome length: ' + str(maxChrLength))

  # Set horizontal scale
  horizontalScale = float(totalWidth) / float(maxChrLength)

  # Draw overall figure with no rectangles
  fig = figure(num=None, figsize=(totalWidthInches + float(totalPaddingPixels + textWidthPixels) / pixelsPerInch + 1, totalHeightInches), dpi=pixelsPerInch, facecolor='w', edgecolor='k')

  colors = ['red', 'orange', 'yellow', 'green', 'blue', 'purple']

  startHeight = spaceHeight
  for i in range(0, len(lines), 2):
    if len(lines[i]) == 0:
      continue
    endHeight = startHeight + chrHeight   
    next = lines[i+1].split()
    nextInt = [int(val) for val in next]
    totalLength = max(nextInt)
    print(totalLength)
    rect = patches.Rectangle((textWidthPixels,totalHeight - endHeight), int(horizontalScale * totalLength), chrHeight,linewidth=1,edgecolor='black',facecolor='none', zorder = 1000, figure = fig)
    fig.patches.extend([rect])
    fig.text(0, float(totalHeight - endHeight - chrHeight) / totalHeight, lines[i], fontsize = 19)
    for j in range(0, len(next), 2):
      start = int(next[j])
      end = int(next[j+1])
      if end - start < 1:
        continue
      rect = patches.Rectangle((textWidthPixels + int(start * horizontalScale),totalHeight - endHeight + 1), int(horizontalScale * (end - start)), chrHeight - 2,linewidth=1,edgecolor=colors[(j/2)%6],facecolor=colors[(j/2)%6], fill = True, zorder = 100, figure = fig, alpha = 0.2)
      fig.patches.extend([rect])
      #print(str(start)+' '+str(end))
    
    startHeight += chrHeight + spaceHeight
  plt.savefig(ofn)
