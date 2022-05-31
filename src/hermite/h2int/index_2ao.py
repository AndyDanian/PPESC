n = 4

i = 0
j = 0

stop_out = False
# while stop_out == False:

#     stop_into = False
#     k = 0
#     l = 0
#     while stop_into == False:
        
#         print(f"({i},{j}|{k},{l})")
         
#         if j > l and k > l:
#             l += 1
#         elif k > j:
#             j += 1
#             stop_into = True
#         elif i == j == k == l:
#             i += 1
#             j = 0
#             stop_into = True

#         if i > k:
#             k += 1        

#         if i == j == k == l == 2:
#             print(f"({i},{j}|{k},{l})")
#             stop_out = stop_into = True

# nocc = 1
# nprim = 2
# icount = 0
# i = 1
# j = 1
# a = nocc + 1
# b = nocc + 1
# if nocc == 1:
#     a = 2
# if nocc == 1:
#     b = 2
# a1 = 1
# b1 = 1

# while icount < 1:
#     print(a1,j,b1,i)

#     i = i + 1
#     if (i > nocc):
#         j = j + 1
#         i = 1

#     if (j > nocc):
#         if (a == b):
#             a = a + 1
#             b = nocc + 1

#             a1 = a1 + 1
#             b1 = 1
#         else:
#             b = b + 1
#             b1 = b1 + 1

#         i = 1
#         j = 1

#     if a > nprim: 
#         icount = 2

count = 1
for a in range(4):
    for j in range(a, 4):
        for b in range(a, 4):
            for i in range(a, 4):
                
                print(f"{count} ({a}{j}|{b}{i})")
                count += 1
