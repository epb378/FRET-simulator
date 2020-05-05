from flask import Flask, render_template, request, jsonify
import os
import tempfile
import matplotlib.pyplot as plt
import nedfunctions as nf
UPLOAD_FOLDER = '/home/littleneddyb/fretmeister/static'
app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
from flask import Flask, render_template, request
def remove_prefix(text, prefix):
    if text.startswith(prefix):
        return text[len(prefix):]
    return text
app = Flask(__name__)
a=tempfile.NamedTemporaryFile()
filename=a.name+'.png'
@app.route('/', methods =['GET', 'POST'])
def upload_file(percent='n/a', spectrum='/static/test.png', filename=filename):
   if request.method == 'POST':
      number1 = request.values['number1']
      number2 = request.values['number2']    
      number3 = request.values['number3']
      plt.clf()
      plt.cla()
      plt.close()
      outplqe, image=nf.FRETsimulator(float(number1),float(number2),float(number3))
      plt.hist(image, bins = 35)
      plt.xlim(350,750)
      plt.xlabel('Emission Wavelength (nm)')
      plt.ylabel('Photons Emitted')
      a=tempfile.NamedTemporaryFile()
      filename=remove_prefix(a.name, '/tmp/')
      filename='./static/'+filename+'.png'
      print(filename)
      plt.savefig(filename)

      percent=str(int(100*outplqe))

   return render_template('index.html', spectrum=filename, percent=percent)

#@app.route('/uploader', methods = ['GET', 'POST'])
#def upload_file2():
#   if request.method == 'POST':
#      f = request.files['file']
#      f.save(secure_filename(f.filename))
#      return 'file uploaded successfully'
		


if __name__ == '__main__':
    # This is used when running locally only. When deploying to Google App
    # Engine, a webserver process such as Gunicorn will serve the app. This
    # can be configured by adding an `entrypoint` to app.yaml.
    # Flask's development server will automatically serve static files in
    # the "static" directory. See:
    # http://flask.pocoo.org/docs/1.0/quickstart/#static-files. Once deployed,
    # App Engine itself will serve those files as configured in app.yaml.
    app.run(host='127.0.0.1', port=8080, debug=True)

