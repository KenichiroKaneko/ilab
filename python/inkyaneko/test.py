pos = r"python/inkyaneko/CCS_FLXLP_sensor_position.txt"

with open(pos, mode='r', encoding='utf-8') as f:
    # f.read()
    print(f.readline())
    l_strip = [s.strip() for s in f.readlines()]

r = []
z = []

for row in l_strip:
    list = row.split('\t')

    r.append(float(list[0]))
    z.append(float(list[1]))

print(r, z)
