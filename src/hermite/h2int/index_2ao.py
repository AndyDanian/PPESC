n = 4

i = 0
j = 0

stop_out = False
while stop_out == False:

    stop_into = False
    k = 0
    l = 0
    while stop_into == False:
        
        print(f"({i},{j}|{k},{l})")
         
        if j > l and k > l:
            l += 1
        elif k > j:
            j += 1
            stop_into = True
        elif i == j == k == l:
            i += 1
            j = 0
            stop_into = True

        if i > k:
            k += 1        

        if i == j == k == l == 2:
            print(f"({i},{j}|{k},{l})")
            stop_out = stop_into = True

# i = 1
# j = 1
# a = 1
# b = 1
      
# icount = 1
# irow = 0
# iprint = 0

# while irow < 4:
#     irow = irow + 1

#     icol = 0

#     if i > 4:
#         j = j + 1
#         i = 1

#     if j > 4:
#         j = 1

#     while icol < 4:
#         icol = icol + 1
#         #Contrucción indices     

#         print(icount,a,b,j,i)

#         icount = icount + 1

#         i = i + 1
#         if i > 4:
#             j = j + 1
#             i = 1

#         if j > 4:
#             if a == b or b == 4:
#                 a = a + 1
#                 b = 1
#                 icol = 5
#             else:
#                 b = b + 1
#                 icol = 5
#             i = 1
#             j = 1