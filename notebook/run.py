import webbrowser, os
filename = 'test.html'
webbrowser.open('file://' + os.path.realpath(filename))
