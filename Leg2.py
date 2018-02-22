import matplotlib.pyplot as plt


y = [(6.123, 8.234, 10.54),(6.123, 8.234, 10.54),(6.123, 8.234, 10.54)]

x = ['hey','yo','ya']
labels = x

for xe, ye in zip(x, y):
    plt.scatter([xe] * len(ye), ye)


ticks = []
for i in range(len(x)):
    ticks.append(i)


print(ticks)
plt.xticks(ticks, labels, rotation='vertical')
plt.axes().set_xticklabels(x)


plt.show()



def yeye():
    y = [(6.123, 8.234, 10.54),(6.123, 8.234, 10.54),(6.123, 8.234, 10.54)]

    x = [('hey'),('yo'),('ya')]

    print(len(y))
    print(len(x))

    plotmepls(x,y)
    #haha(x1,y,i)

def plotmepls(mainpoint, diffpoints):
    #x = [(mainpoint[i])]
    #y = [(diffpoints[i][0],diffpoints[i][1],diffpoints[i][2])]
    x = mainpoint
    y = diffpoints
    print(len(y))
    print(len(x))

    #fig, ax = plt.subplots()

    #for xe, ye in zip(x, y):
    #    ax.scatter([xe] * len(ye), ye)
    #plt.plot(x, y, s=60, c='red', marker='^')

    #plt.xticks([number], [number])

    #ax.set_ylim(0, 10)
    #plt.xticks([j])
    #plt.axes().set_xticklabels([1])

    #ax.set_xticks(mainpoint)
    plt.scatter(x, y)

    plt.show()

def haha(mainpoint, diffpoints, number):
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker

    x = [(mainpoint,mainpoint,mainpoint)]
    y = [(diffpoints[0], diffpoints[1], diffpoints[2])]

    tick_spacing = 1

    fig, ax = plt.subplots(1,1)
    ax.scatter(x, y)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
    plt.show()