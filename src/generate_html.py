"""
generate_html.py: Generate the HTML file from our plot images.
"""

import logging


class ST_HTML_Generator:
    def __init__(self, para_list, plot_filepaths, static=True):
        self.image_key_list = para_list[0]  # List of statistics variables to look for
        self.header_info = para_list[1]  # The webpage's title
        self.input_para = para_list[2]  # List of the input parameters used for these statistics
        self.static = static  # Static vs. dynamic webpage boolean
        self.plot_filepaths = plot_filepaths
        self.prg_name = self.input_para["prg_name"]  # Program name
        self.html_writer = None

        if len(self.input_para["input_files"]) > 1:
            self.more_input_files = True
        else:
            self.more_input_files = False

    def generate_header(self):
        """Format the header of the HTML file with the title and CSS."""
        html_filepath = self.input_para["output_folder"] + '/' + self.input_para["out_prefix"] + ".html"
        self.html_writer = open(html_filepath, 'w', encoding='utf-8')
        self.html_writer.write("<html>")
        self.html_writer.write("<head>")
        self.html_writer.write("<title>")
        self.html_writer.write(self.header_info)
        self.html_writer.write("</title>")
        self.html_writer.write('''<style  type="text/css">
 @media screen {
  div.summary {
    width: 18em;
    position:fixed;
    top: 3em;
    margin:1em 0 0 1em;
  }
  
  div.main {
    display:block;
    position:absolute;
    overflow:auto;
    height:auto;
    width:auto;
    top:4.5em;
    bottom:2.3em;
    left:18em;
    right:0;
    border-left: 1px solid #CCC;
    padding:0 0 0 1em;
    background-color: white;
    z-index:1;
  }
  
  div.header {
    background-color: #EEE;
    border:0;
    margin:0;
    padding: 0.5em;
    font-size: 200%;
    font-weight: bold;
    position:fixed;
    width:100%;
    top:0;
    left:0;
    z-index:2;
  }

  div.footer {
    background-color: #EEE;
    border:0;
    margin:0;
	padding:0.5em;
    height: 1.3em;
	overflow:hidden;
    font-size: 100%;
    font-weight: bold;
    position:fixed;
    bottom:0;
    width:100%;
    z-index:2;
  }
  
  img.indented {
    margin-left: 3em;
  }
 }
 
 @media print {
	img {
		max-width:100% !important;
		page-break-inside: avoid;
	}
	h2, h3 {
		page-break-after: avoid;
	}
	div.header {
      background-color: #FFF;
    }
	
 }
 
 body {    
  font-family: sans-serif;   
  color: #000;   
  background-color: #FFF;
  border: 0;
  margin: 0;
  padding: 0;
  }
  
  div.header {
  border:0;
  margin:0;
  padding: 0.5em;
  font-size: 200%;
  font-weight: bold;
  width:100%;
  }    
  
  #header_title {
  display:inline-block;
  float:left;
  clear:left;
  }
  #header_filename {
  display:inline-block;
  float:right;
  clear:right;
  font-size: 50%;
  margin-right:2em;
  text-align: right;
  }

  div.header h3 {
  font-size: 50%;
  margin-bottom: 0;
  }
  
  div.summary ul {
  padding-left:0;
  list-style-type:none;
  }
  
  div.summary ul li img {
  margin-bottom:-0.5em;
  margin-top:0.5em;
  }
	  
  div.main {
  background-color: white;
  }
      
  div.module {
  padding-bottom:1.5em;
  padding-top:1.5em;
  }
	  
  div.footer {
  background-color: #EEE;
  border:0;
  margin:0;
  padding: 0.5em;
  font-size: 100%;
  font-weight: bold;
  width:100%;
  }


  a {
  color: #000080;
  }

  a:hover {
  color: #800000;
  }
      
  h2 {
  color: #800000;
  padding-bottom: 0;
  margin-bottom: 0;
  clear:left;
  }

  table { 
  margin-left: 3em;
  text-align: center;
  }
  
  th { 
  text-align: center;
  background-color: #000080;
  color: #FFF;
  padding: 0.4em;
  }      
  
  td { 
  font-family: monospace; 
  text-align: left;
  background-color: #EEE;
  color: #000;
  padding: 0.4em;
  }

  img {
  padding-top: 0;
  margin-top: 0;
  border-top: 0;
  }
  
  p {
  padding-top: 0;
  margin-top: 0;
  }
  
  li {
  margin: 10px 0;
  }
      </style>''')

        self.html_writer.write("</head>")
        self.html_writer.write("<body>")
        self.html_writer.write('''<div class="header"> 
         <div id="header_title">{} Report</div>
         <div id="header_filename">
             <script> document.write(new Date().toLocaleDateString()); </script>
      '''.format(self.prg_name))
        self.html_writer.write('''</div></div>''')

    def generate_left(self):
        """Generate the left section of the HTML file with the links to the
        right section."""
        # Add the summary section of links
        self.html_writer.write('<div class="summary">');
        self.html_writer.write('<h2>Summary</h2>')
        self.html_writer.write('<ul>')

        # Add links to the right sections
        key_index = 0
        for plot_key in self.image_key_list:
            self.html_writer.write('<li>')

            self.html_writer.write(
                '<a href="#lrst' + str(key_index) + '">' + self.plot_filepaths[plot_key]['title'] + '</a>')
            key_index += 1
            self.html_writer.write('</li>')

        # Add the input files section link
        self.html_writer.write('<li>')
        self.html_writer.write('<a href="#lrst' + str(key_index) + '">Input File List</a>')
        key_index += 1
        self.html_writer.write('</li>')
        
        self.html_writer.write("</ul>")
        self.html_writer.write('</div>')

    def generate_right(self):
        """Generate the right section of the HTML file with the plots and tables."""
        self.html_writer.write('<div class="main">')
        key_index = 0
        for plot_key in self.image_key_list:
            self.html_writer.write('<div class="module">')
            self.html_writer.write(
                '<h2 id="lrst' + str(key_index) + '">' + self.plot_filepaths[plot_key]['description'] + '</h2><p>')

            # Add the figures
            if plot_key == "basic_st" or plot_key == "base_mods":
                # Add the HTML tables
                self.html_writer.write(self.plot_filepaths[plot_key]['detail'])

            else:
                # Add the dynamic plots
                try:
                    dynamic_plot = self.plot_filepaths[plot_key]['dynamic']
                    self.html_writer.write(dynamic_plot)

                except KeyError:
                    logging.error("Missing dynamic plot for %s", plot_key)

            self.html_writer.write('</div>')

            key_index += 1

        self.html_writer.write('<div class="module">')
        self.html_writer.write('<h2 id="lrst' + str(key_index) + '">File count = ' + str(
            len(self.input_para["input_files"])) + '</h2><p>')
        for _af in self.input_para["input_files"]:
            self.html_writer.write("<br/>" + _af)
        self.html_writer.write('</p></div>')
        key_index += 1

        self.html_writer.write('</div>')
    
    def generate_left_signal_data(self, read_names):
        """Generate the left section of the HTML file with the links to the right section."""
        self.html_writer.write('<div class="summary">');
        self.html_writer.write('<h2>Summary</h2>')
        self.html_writer.write('<ul>')

        # Add the summary table section link
        url_index = 0
        self.html_writer.write('<li>')
        self.html_writer.write(
            '<a href="#lrst' + str(url_index) + '">' + self.plot_filepaths["basic_st"]['title'] + '</a>')
        url_index += 1

        # Add the read plot section list
        for read_name in read_names:
            self.html_writer.write('<li>')
            self.html_writer.write(
                '<a href="#lrst' + str(url_index) + '">' + read_name + '</a>')
            url_index += 1

        # Add the input files section link
        self.html_writer.write('<li>')
        self.html_writer.write('<a href="#lrst' + str(url_index) + '">Input Files</a>')
        url_index += 1
        self.html_writer.write('</li>')
        self.html_writer.write("</ul>")
        self.html_writer.write('</div>')

    # Generate tables and plots in the right section
    def generate_right_signal_data(self, read_names, signal_plot):

        self.html_writer.write('<div class="main">')

        # Add the summary table section
        url_index = 0
        self.html_writer.write('<div class="module">')
        self.html_writer.write(
            '<h2 id="lrst' + str(url_index) + '">' + self.plot_filepaths["basic_st"]['description'] + '</h2><p>')
        self.html_writer.write(self.plot_filepaths["basic_st"]['detail'])
        url_index += 1

        # Add the read plot section
        for read_name in read_names:
            self.html_writer.write('<div class="module">')

            # Set the description
            description_text = "ONT Basecall Signal"
            self.html_writer.write(
                '<h2 id="lrst' + str(url_index) + '">' + description_text + '</h2><p>')

            # Add the plot
            html_plot = signal_plot.get(read_name)
            self.html_writer.write(html_plot)
            self.html_writer.write('</div>')
            url_index += 1

        # Add the input files section
        self.html_writer.write('<div class="module">')
        self.html_writer.write('<h2 id="lrst' + str(url_index) + '">File count = ' + str(
            len(self.input_para["input_files"])) + '</h2><p>')
        for input_file_str in self.input_para["input_files"]:
            self.html_writer.write("<br/>" + input_file_str)
        self.html_writer.write('</p></div>')
        self.html_writer.write('</div>')

    def generate_end(self):
        self.html_writer.write(
            '<div class="footer"> Generated by <a href="https://github.com/WGLab/{}">{}</a> </div>'.format(
                self.prg_name, self.prg_name))

        self.html_writer.write("</body>")
        self.html_writer.write("</html>")
        self.html_writer.close()

    # Main function for generating the HTML.
    def generate_st_html(self, signal_plots=False):
        if signal_plots:
            self.generate_header()
            # Get the signal plots
            signal_plots = self.plot_filepaths["ont_signal"]['dynamic']
            read_names = signal_plots.keys()
            self.generate_left_signal_data(read_names)
            self.generate_right_signal_data(read_names, signal_plots)
            self.generate_end()
        else:
            # Format base QC
            self.generate_header()
            self.generate_left()
            self.generate_right()
            self.generate_end()


