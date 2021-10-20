import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches

class colorDf(pd.DataFrame):
    def __init__(self, df, colorDict):
        super().__init__(df)
        self.colorDict = colorDict

def multiLandscape(dfs, figsize = (220,50)):
    # calculate amount of rows for all dataframe
    allRow = 0
    maxNumColor = 1
    for df in dfs: 
        allRow += len(df.index)
        maxNumColor = max(5*len(df.colorDict.keys()),maxNumColor)

    plt.subplots(figsize = figsize)
    grid = plt.GridSpec(allRow, int(len(dfs[0].columns)+maxNumColor)+1)

    # constants for rectangle 
    xSepDist = 0.5
    ySepDist = 1
    width = 2
    height = 10
    
    rowStart = 0

    for df in dfs:
        # new subplot for current df
        rowEnd = rowStart + len(df.index)
        ax = plt.subplot(grid[rowStart:rowEnd,0:len(df.columns)])
        ax.axis([0,2.5*(len(df.columns)),0,11*len(df.index)])

        xStartPoint = 0
        yStartPoint = 0
        yPoint = yStartPoint

        rIndex = df.index.to_list()
        rIndex.reverse()
        for i in rIndex: #each row (gene)
            xPoint = xStartPoint
            for j in df.columns: #draw one row #number of columns
                tempRectangle = patches.Rectangle((xPoint,yPoint),width,height,linewidth=3,edgecolor='black',facecolor=df.colorDict[str(df.loc[i,j])])
                xPoint = xPoint + width + xSepDist
                ax.add_patch(tempRectangle)
            yPoint = yPoint + height + ySepDist

        ax.set_frame_on(False) #remove the lines
        plt.tick_params(top=False, bottom=False, left=False, right=False, labelleft=True, labelbottom=True) #remove ticks
        xGap = (xPoint - xStartPoint)/len(df.columns)
        yGap = (yPoint - yStartPoint)/len(df.index)

        ax.set_yticks([yStartPoint+(0.5*height)+ySepDist+(y*yGap) for y in range(len(df.index))])
        subIndex = df.index.tolist()
        subIndex.reverse()
#         ax.set_yticklabels(subIndex) # draw yticklabels
#         yPoint = yPoint + height
        ax.set_xticks([xStartPoint+(0.5*width)+xSepDist+(x*xGap) for x in range(len(df.columns))])
        ax.set_xticklabels(['' for x in range(len(df.columns))]) # hide xticklabels

#         ax2 = plt.subplot(grid[rowStart:rowEnd,len(df.columns)+1:len(df.columns)+1+maxNumColor])
#         ax2.axis([0,2.5*(maxNumColor),0,11*len(df.index)])
#         xPoint = xStartPoint
#         for i in df.colorDict.keys(): #draw one row #number of columns
#             tempRectangle = patches.Rectangle((xPoint,yPoint-height-ySepDist),width,height,linewidth=0.5,edgecolor='black',facecolor=df.colorDict[i])
#             xPoint = xPoint + 5*width + 5*xSepDist
#             ax2.add_patch(tempRectangle)
            
#         ax2.set_xticklabels(['' for x in range(len(df.columns))]) # hide xticklabels
#         ax2.set_yticklabels(['' for y in range(len(df.columns))]) # hide yticklabels   
#         ax2.set_frame_on(False) #remove the lines
        
        rowStart = rowEnd
    
    ax.set_xticklabels(dfs[0].columns.tolist(), rotation = 90) # draw xticklabels for the last df