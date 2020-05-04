
# coding: utf-8

# In[1]:


import numpy as np 
import pandas as pd


# In[2]:


def get_blocks(file_name):
    '''open a nef file and get the blocks from the file.
    Input: a string as file name.
    Retrun: a list of all the blocks(also a list). 
    In one block list, lines are elements in the list.
    '''
    f = open(file_name, 'r')
    lines = list(f)
    all_blocks = []
    block = []
    for i in range(len(lines)):
        if lines[i].startswith('   save_nef_'):
          #  block = []
            block.append(lines[i])
            for j in range(i + 1, len(lines)):
                block.append(lines[j])
                if lines[j].startswith('   save_\n'):
                    all_blocks.append(block)
                    block = []
                    break   
    return all_blocks


# In[3]:


all_blocks = get_blocks('CCPN_Commented_Example.nef')


# In[4]:


#use dictionary to store all the blocks
index_all = []
item_all = []

for block in all_blocks:
    #use dictionary to store all the info in one block
    index = []
    item = []
    #first, basic info of the block
    basic_info_1 = []
    basic_info_2 = []
    for i in range(len(block)):
        if block[i].startswith('      _nef_'):
            line = block[i].split()
            name = line[0].split('.')[0]
            basic_info_1.append(line[0])
            basic_info_2.append(line[1])
    df = pd.DataFrame(np.array([basic_info_2]), columns=basic_info_1 )
    index.append(name)
    item.append(df)

    #inside the block, what remains is loops. To get the loop:
    all_loops = []
    loop = []
    for i in range(len(block)):
        if block[i].startswith('      loop_\n'):
            for j in range(i + 1, len(block)):
                loop.append(block[j])
                if block[j].startswith('      stop_\n'):
                    loop = loop[:len(loop) - 1]
                    all_loops.append(loop)
                    loop = []
                    break   

    for i in range(len(all_loops)):
        #choose one loop
        loop = all_loops[i]
        #get the column name list of this loop
        col = []
        for line in loop:
            line = line.strip()
            col.append(line)
            if not line:
                col = col[:len(col) - 1]
                break
        #get the data in this loop
        for i in range(len(loop)):
            if loop[i].startswith('\n'):
                data = []
                for j in range(i + 1, len(loop)):
                    data.append(loop[j])
        #make a 2d empty list to store the data
        col_num = len(data[0].split())
        empty_lists = [ [] for i in range(col_num) ]
        for line in data:
            col_num = len(line.split())
            row = line.split()
            for i in range(len(empty_lists)):
                li = empty_lists[i]
                element = row[i]
                li.append(element)
        #make a dataframe to store the data
        loop_data = pd.DataFrame(np.array(empty_lists).T, columns=col) 
        loop_name = col[0].split('.')[0]
        index.append(loop_name)
        item.append(loop_data)
    block_data = dict(zip(index,item))
    block_name = block[0].split()[0]
    index_all.append(block_name)
    item_all.append(block_data)

nef = dict(zip(index_all,item_all))
nef.keys()


# In[7]:


nef['save_nef_nmr_meta_data'].keys()


# In[8]:


nef['save_nef_nmr_meta_data']['_nef_program_script']


# In[ ]:


#Ends here.


# In[ ]:


txt = "hello, my name is Peter, I am 26 years old"

x = txt.split(", ")

print(x)


# In[ ]:


pp.to_string()


# In[ ]:


D = {'a': pp}
#dictionry would be good to map thge blocks
#inside the blocks dictionary would be good to map loops
D


# In[ ]:


import io 
data_string = """Letters, Numbers
                 a, 1
                 b, 2
                 c, 3"""

data = io.StringIO(data_string)
df = pd.read_csv(data, sep=",")
df


# In[ ]:


#this is the method to get blocks
all_blocks = []
for i in range(len(lines)):
    if lines[i].startswith('   save_nef_'):
        block = []
        block.append(lines[i])
        for j in range(i +1 , len(lines)):
            block.append(lines[j])
            if lines[j].startswith('   save_\n'):
                all_blocks.append(block)
            #    block = []
                break
                
                
                
                


# In[ ]:


with open("demofile2.txt", "w") as fo:
    for i in all_blocks[0]:
        fo.write(i)


# In[ ]:


i = 1
outfile = open(f'outfile{i}.txt', 'w')
for line in infile:
    outfile.write(line.strip())
    if is_last_line(line):
        i += 1
        outfile = open(f'outfile{i}.txt', 'w')
outfile.close()


# In[ ]:


df2 = pd.DataFrame(np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]), columns=['a', 'b', 'c'])
df2


# In[ ]:


line = ''
if line:
    print('y')


# In[ ]:


fo = open('CCPN_Commented_Example.nef', 'r')


# In[ ]:


a  = fo.readlines()


# In[ ]:


a


# In[ ]:


whole_file = []
for line in a:
    line = line.strip()
    whole_file.append(line)


# In[ ]:


for line in whole_file:
    if line.startswith('save_'):
        print(line)

