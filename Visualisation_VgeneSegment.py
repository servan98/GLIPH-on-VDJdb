# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 15:37:09 2020

@author: S.J. van den Brink

In order to run the visualisation for your own cluster, first run the R file and the gliph clustering.
Then change the files to the correct files in line 22,28,33,92,118 and 132.
To move to a different cluster, only line 32 and 125 need to be changed
"""

#import packages
from pyvis.network import Network
import pandas as pd

#create empty network
klr_net = Network(height="75%", width="100%", bgcolor="#ffffff", font_color="#004a47")

# set the physics layout of the network
klr_net.barnes_hut()

klr = pd.read_table("D:/gliph-1.0/gliph/bin/DATAz-clone-network.txt", header = None)
sources = klr[0]
targets = klr[1]
edge_data = zip(sources, targets)

#open HLA data made in R
with open("D:/gliph-1.0/gliph/bin/VgeneUsage_DATA.txt") as b:
  Vs = b.readlines()
b.close()

#open cluster data made in R and make a list out of it
with open("D:/gliph-1.0/gliph/bin/cluster457.txt") as f:
  cl18 = f.readline()
f.close()
cl18 = cl18.split()

x=0
for get in Vs:
    x+=1
print(x)


for e in edge_data:     #for every source and target after gliph
    src = e[0]  #assign source to src
    dst = e[1]  #assign target/destination to dst
    if src in cl18 and dst in cl18: #check if source and data are both present in convergencegroup (18)
        for i in range(0,x):   #itterate over 0-35 for every HLA in HLAs
            if src in Vs[i]:  #if source in one of the HLA lists..
                klr_net.add_node(src, src, title=src)   #add source node to prob_net with color assigned to i
        for i in range(0,x):   #itterate over 0-35 for every HLA in HLAs
            if dst in Vs[i]:  #if target in one of the HLA lists..
                klr_net.add_node(dst, dst, title=dst) #add target node to prob_net with color assigned to i     
        klr_net.add_edge(src, dst) #add edge between source and target
neighbor_map = klr_net.get_adj_list() #create map with neighbours

# add neighbor data to node hover data
for node in klr_net.nodes:
    node["title"] += " Neighbors:<br>" + "<br>".join(neighbor_map[node["id"]])
    node["value"] = len(neighbor_map[node["id"]])

#Make a dictionary containing the occurence of every HLA-type
occ = {}
for point in klr_net.nodes: #for every node
    for i in range(0,x): # create a loop over 35 numbers(enough to cover all HLAs in this case)
        if point["id"] in Vs[i]: #if node sequence in the list of HLA[i]
            name = Vs[i].split()[0] #assign the HLA-type to name
            if name in occ.keys(): #if name is present in occ as key
                occ[name] = occ[name]+1 #add 1 to the value of this key
            else: #otherwise, create key with value 1
                occ[name] = 1

#list of differentiatable colours
#           rood        green       blauw       geel        paars       lichtblauw      roze        oranje     zwart    lichtpaars   donkergroen
Hcolor = {0:"#ff0000",1:"#22ff00",2:"#002fff",3:"#f6ff00",4:"#9700bd",5:"#00fffb",6:"#ff00e6",7:"#ffae00",8:"#000000",9:"#d191ff",10:"#228713",}
#assign a readable colour to each HLA type from the list above if occurence is greater than 2, else make it grey
Ncolors = {}
occ = {k: v for k, v in sorted(occ.items(), key=lambda item: item[1], reverse=True)}
z = 0
for HLAn, amm in occ.items():
    if amm > 4:
        Ncolors[HLAn] = Hcolor[z]
        z += 1
    else:
        Ncolors[HLAn] = "#bababa"

#make an empty network
prob_net = Network(height="75%", width="100%", bgcolor="#ffffff", font_color="#004a47")

# set the physics layout of the network
prob_net.barnes_hut()
prob = pd.read_table("D:/gliph-1.0/gliph/bin/DATAz-clone-network.txt", header = None)
sources = prob[0]
targets = prob[1]
edge_data = zip(sources, targets)

for e in edge_data:     #for every source and target after gliph
    src = e[0]  #assign source to src
    dst = e[1]  #assign target/destination to dst
    if src in cl18 and dst in cl18: #check if source and data are both present in convergencegroup (18)
        for i in range(0,x):   #itterate over 0-35 for every HLA in HLAs
            if src in Vs[i]:  #if source in one of the HLA lists..
                name = Vs[i].split()[0] #get the HLA-type of this source
                prob_net.add_node(src, src, title=src, color=Ncolors[name])   #add source node to prob_net with color assigned to the hla-type in Ncolors
        for i in range(0,x):   #itterate over 0-35 for every HLA in HLAs
            if dst in Vs[i]:  #if target in one of the HLA lists..
                name = Vs[i].split()[0] #get the HLA-type of this target
                prob_net.add_node(dst, dst, title=dst, color=Ncolors[name]) #add target node to prob_net with color assigned to the hla-type in Ncolors    
        prob_net.add_edge(src, dst) #add edge between source and target
neighbor_map = prob_net.get_adj_list()

# add neighbor data to node hover data
for node in prob_net.nodes:
    node["title"] += " Neighbors:<br>" + "<br>".join(neighbor_map[node["id"]])
    node["value"] = len(neighbor_map[node["id"]])

#save graph as TCR.html
prob_net.save_graph("D:/gliph-1.0/gliph/bin/TCR.html")

#make a legend in html language
legenda = "<div> \n<h3>Legend HLA types</h3> \n V gene usage : occurence"
for key, value in Ncolors.items():
    legenda += "<p><font color=" + value + ">" + key + "=" + str(occ[key]) + "</p> \n"
legenda+="</div> \n"

#open the previously saved TCR.html, remove the "</body> </html>" add the legend and add "</body> </html>" again. then write to the appropriate file 
with open('D:/gliph-1.0/gliph/bin/TCR.html', 'r') as f:
    html = f.read()
htmls = html[-0:-15]
htmll = htmls+legenda
htmlf = htmll + html[-15:]
H= open("D:/gliph-1.0/gliph/bin/TCRL457V.html","w+")
H.write(htmlf)
H.close()
#list of clusters i have used in visualisation: 18,2,139,111,33,113,591,457,1370,2534
