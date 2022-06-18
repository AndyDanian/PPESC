n = 44 - 1

stop = False
i = j = a = b = 0
# i j a b
# i j b a
# a b i j

count = 0
while not stop:

    if a > i:
        i += 1
        j = a = b = 0

    #print(f"ijkl {i + 1}  {j + 1}  {a + 1}  {b + 1}")
    if i == n  and j == n  and a == n  and b == n :
        stop = True

    if b < a:
        b += 1
    elif i == a and a > j:
        j += 1
        a = b = 0
    else:
        a += 1
        b = 0
    count += 1

print("# : ",count)
    # if a > j and b + 1 >= i:
    #     j += 1
    #     b = 0
    #     if a - j == 0:
    #         a = 0
    # elif b + 1 < i and b < a:
    #     b += 1
    # elif j == a and a > b:
    #     b += 1
    # elif b > j:
    #     j += 1
    #     a = b = 0
    # else:
    #     a += 1
    #     if j == b and b < a:
    #         j = b = 0
    #     elif b < a:
    #         b = 0