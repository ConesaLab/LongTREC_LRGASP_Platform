from app import app, db
from flask import render_template, url_for, flash, redirect, request, render_template, request, send_file, Flask, render_template, jsonify, Response, send_from_directory
from forms import RegistrationForm, LoginForm, SubmissionForm
from models import User, Submission
from flask_login import login_user, current_user, logout_user, login_required
import os, shutil, time
from werkzeug.utils import secure_filename
from utils import allowed_file, evaluate_submission
import subprocess
import glob
from threading import Thread
import json
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

logger.info("This is an info message")
logger.error("This is an error message")

print("Hello world", flush=True)
print("Hello world", flush=True)

progress = 0
status_message = "Booting Script..."
terminal_output = []

@app.route("/")
@app.route("/home")
def home():
    return render_template('index.html')

@app.route("/register", methods=['GET', 'POST'])
def register():
    if current_user.is_authenticated:
        return redirect(url_for('dashboard'))
    form = RegistrationForm()
    if form.validate_on_submit():
        # Hash the password here in production!
        user = User(username=form.username.data, email=form.email.data, password=form.password.data)
        db.session.add(user)
        db.session.commit()
        flash('Account created! You can now log in.', 'success')
        return redirect(url_for('login'))
    return render_template('register.html', title='Register', form=form)

@app.route("/login", methods=['GET', 'POST'])
def login():
    if current_user.is_authenticated:
        return redirect(url_for('dashboard'))
    form = LoginForm()
    if form.validate_on_submit():
        user = User.query.filter_by(email=form.email.data).first()
        if user and user.password == form.password.data:
            # For security, implement password hashing!
            login_user(user, remember=form.remember.data)
            next_page = request.args.get('next')
            return redirect(next_page) if next_page else redirect(url_for('dashboard'))
        else:
            flash('Login unsuccessful. Please check email and password', 'danger')
    return render_template('login.html', title='Login', form=form)

@app.route("/logout")
def logout():
    logout_user()
    return redirect(url_for('home'))

@app.route("/dashboard")
@login_required
def dashboard():
    submissions = Submission.query.filter_by(author=current_user)
    return render_template('dashboard.html', title='Dashboard', submissions=submissions)

@app.route("/submit", methods=['GET', 'POST'])
@login_required
def submit():
    form = SubmissionForm()
    if form.validate_on_submit():
        if form.data_file.data and allowed_file(form.data_file.data.filename):
            filename = secure_filename(form.data_file.data.filename)
            file_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
            form.data_file.data.save(file_path)
            # Evaluate the submission
            evaluation_result = evaluate_submission(file_path)
            # Save submission to the database
            submission = Submission(data_file=filename, evaluation_result=evaluation_result, author=current_user)
            db.session.add(submission)
            db.session.commit()
            flash('Submission successful!', 'success')
            return redirect(url_for('dashboard'))
    return render_template('upload.html', title='Submit Prediction', form=form)


@app.route("/submission/<int:submission_id>")
@login_required
def submission(submission_id):
    submission = Submission.query.get_or_404(submission_id)
    if submission.author != current_user:
        abort(403)
    return render_template('submission.html', title='Submission Detail', submission=submission)

@app.route("/data")
def data():
    return render_template('data.html', title='Data')

@app.route("/challenge1", methods=['GET', 'POST'])
def challenge1():
    return render_template("challenge1_final.html")

@app.route("/run_script_challenge1", methods=['POST'])
def run_script_challenge1():
    print('dikke bictor')
    progress = 0
    status_message = "Starting Challenge 1..."
    terminal_output = []

    def replace_directory(directory_path):
        if os.path.exists(directory_path):
            shutil.rmtree(directory_path)
        os.makedirs(directory_path)

    replace_directory("uploads")
    replace_directory("uploads/transcriptome_file1")
    replace_directory("uploads/transcriptome_file2")
    replace_directory("sqanti_results/results_file1")
    replace_directory("sqanti_results/results_file2")
    replace_directory("uploads/coverage_files")
    replace_directory("uploads/coverage_files2")

    organism = request.form.get('organism')

    # Get file and form data from the request
    transcriptome_file = request.files.get('file')
    transcriptome_path = os.path.join(app.config['UPLOAD_FOLDER'], 'transcriptome_file1', transcriptome_file.filename)
    transcriptome_file.save(transcriptome_path)

    # Get file and form data from the request
    transcriptome_file2 = request.files.get('file-2')
    print('tr_file2:', transcriptome_file2.filename)
    if transcriptome_file2.filename != '':
        print('gekke shit')
        transcriptome_path2 = os.path.join(app.config['UPLOAD_FOLDER'], 'transcriptome_file2', transcriptome_file2.filename)
        transcriptome_file2.save(transcriptome_path2)
    else:
        print('tr_file is NA')
        transcriptome_path2 = 'NA'

    annotation_file = request.form.get('annotation')
    if annotation_file == 'default':
        annotation_path = 'LRGASP_DATA'
    else:
        annotation_file = request.files.get('annotation_file')
        annotation_path = os.path.join(app.config['UPLOAD_FOLDER'], 'transcriptome_file1', annotation_file.filename)
        annotation_file.save(annotation_path)

    annotation_file2 = request.form.get('annotation2')
    print('anno2:', annotation_file2)
    if annotation_file2 != None:
        if annotation_file2 == 'default':
            annotation_path2 = 'LRGASP_DATA'
        else:
            annotation_file2 = request.files.get('annotation_file-2')
            annotation_path2 = os.path.join(app.config['UPLOAD_FOLDER'], 'transcriptome_file2', annotation_file2.filename)
            annotation_file2.save(annotation_path2)
    else:
        annotation_path2 = 'NA'

    reference_file = request.files.get('reference_file')
    if reference_file and reference_file.filename != '':
        reference_path = os.path.join(app.config['UPLOAD_FOLDER'], 'transcriptome_file1', reference_file.filename)
        reference_file.save(reference_path)
    else:
        reference_path = 'LRGASP_DATA'

    reference_file2 = request.files.get('reference_file-2')
    if reference_file2 != None:
        if reference_file2 and reference_file2.filename != '':
            reference_path2 = os.path.join(app.config['UPLOAD_FOLDER'], 'transcriptome_file2', reference_file2.filename)
            reference_file2.save(reference_path2)
        else:
            reference_path2 = 'LRGASP_DATA'
    else: reference_path2 = 'NA'

    coverage = request.form.get('coverage_directory')
    if coverage == 'custom':
        coverage_files = request.files.getlist('coverage_directory[]')
        coverage_dir = coverage_files[0].filename.split('/')[0]
        replace_directory(os.path.join("uploads/coverage_files", coverage_dir))
        for cov_file in coverage_files:
            coverage_path = os.path.join(app.config['UPLOAD_FOLDER'], 'coverage_files', cov_file.filename)
            cov_file.save(coverage_path)
    else:
        coverage_dir = 'NA'

    coverage2 = request.form.get('coverage_directory-2')
    if coverage2 == 'custom':
        coverage_files2 = request.files.getlist('coverage_directory_2[]')
        coverage_dir2 = coverage_files2[0].filename.split('/')[0]
        replace_directory(os.path.join("uploads/coverage_files2", coverage_dir2))
        for cov_file in coverage_files2:
            coverage_path2 = os.path.join(app.config['UPLOAD_FOLDER'], 'coverage_files2', cov_file.filename)
            cov_file.save(coverage_path2)
    else:
        coverage_dir2 = 'NA'
        coverage2 = 'NA'

    cage_file = request.form.get('cage')
    print('cage_file:', cage_file)
    if cage_file == 'default':
        cage_path = 'LRGASP_DATA'
    elif cage_file == 'reference':
        cage_path = 'reference'
    else:
        cage_file = request.files.get('cage_file')
        cage_path = os.path.join(app.config['UPLOAD_FOLDER'], 'transcriptome_file1', cage_file.filename)
        cage_file.save(cage_path)

    cage_file2 = request.form.get('cage2')
    print('cage_file2:', cage_file2)
    if cage_file2 != None:
        if cage_file2 == 'default':
            cage_path2 = 'LRGASP_DATA'
        elif cage_file2 == 'reference':
            cage_path2 = 'reference'
        else:
            cage_file2 = request.files.get('cage_file2')
            cage_path2 = os.path.join(app.config['UPLOAD_FOLDER'], 'transcriptome_file2', cage_file2.filename)
            cage_file2.save(cage_path2)
    else:
        cage_path2 = 'NA'

    quant_file = request.form.get('quant')
    if quant_file == 'default':
        quant_path = 'LRGASP_DATA'
    elif quant_file == 'reference':
        quant_path = 'reference'
    else:
        quant_file = request.files.get('quant_file')
        quant_path = os.path.join(app.config['UPLOAD_FOLDER'], 'transcriptome_file1', quant_file.filename)
        quant_file.save(quant_path)

    quant_file2 = request.form.get('quant2')
    print('quant_file2', quant_file2)
    if quant_file2 != None:
        if quant_file2 == 'default':
            quant_path2 = 'LRGASP_DATA'
        elif quant_file2 == 'reference':
            quant_path2 = 'reference'
        else:
            quant_file2 = request.files.get('quant_file2')
            quant_path2 = os.path.join(app.config['UPLOAD_FOLDER'], 'transcriptome_file2', quant_file2.filename)
            quant_file2.save(quant_path2)
    else:
        quant_path2 = 'NA'

    poly_A_file = request.form.get('polyA')
    if poly_A_file == 'default':
        poly_A_path = 'LRGASP_DATA'
    elif poly_A_file == 'NA':
        poly_A_path = 'NA'
    else:
        poly_A_file = request.files.get('polyA_file')
        poly_A_path = os.path.join(app.config['UPLOAD_FOLDER'], 'transcriptome_file1', poly_A_file.filename)
        poly_A_file.save(poly_A_path)

    poly_A_file2 = request.form.get('polyA2')
    print('poly_A_file2', poly_A_file2)
    if poly_A_file2 != None:
        if poly_A_file2 == 'default':
            poly_A_path2 = 'LRGASP_DATA'
        elif poly_A_file2 == 'NA':
            poly_A_path2 = 'NA'
        else:
            poly_A_file2 = request.files.get('polyA_file2')
            poly_A_path2 = os.path.join(app.config['UPLOAD_FOLDER'], 'transcriptome_file2', poly_A_file2.filename)
            poly_A_file2.save(poly_A_path2)
    else:
        poly_A_path2 = 'NA'


    tool = request.form.get('tool')
    platform = request.form.get('platform')
    platform2 = request.form.get('platform-2')
    if platform2 == None:
        platform2 = 'NA'
    library_preparation = request.form.get('library_preparation')
    library_preparation2 = request.form.get('library_preparation-2')
    if library_preparation2 == None:
        library_preparation2 = 'NA'
    data_category = request.form.get('data_category')
    data_category2 = request.form.get('data_category-2')
    if data_category2 == None:
        data_category2 = 'NA'

    sirv_list = request.files.get('sirv_file')
    if sirv_list:
        # Save the file if it exists
        file_path = os.path.join(app.config['UPLOAD_FOLDER'], 'sirv_list', sirv_list.filename)
        os.makedirs(os.path.dirname(file_path), exist_ok=True)  # Create directory if not exists
        sirv_list.save(file_path)
        sirv_list = os.path.join(app.config['UPLOAD_FOLDER'], 'sirv_list', sirv_list.filename)
    else:
        # If no file is provided, set sirv_list to NA
        sirv_list = "NA"

    sirv_list2 = request.files.get('sirv_file2')
    print('sirv_list2', sirv_list2)
    if sirv_list2:
        # Save the file if it exists
        file_path = os.path.join(app.config['UPLOAD_FOLDER'], 'sirv_list2', sirv_list2.filename)
        os.makedirs(os.path.dirname(file_path), exist_ok=True)  # Create directory if not exists
        sirv_list2.save(file_path)
        sirv_list2 = os.path.join(app.config['UPLOAD_FOLDER'], 'sirv_list2', sirv_list2.filename)
    else:
        # If no file is provided, set sirv_list to NA
        sirv_list2 = "NA"

    ercc_list = request.files.get('ercc_file')
    if ercc_list:
        # Save the file if it exists
        file_path = os.path.join(app.config['UPLOAD_FOLDER'], 'ercc_list', ercc_list.filename)
        os.makedirs(os.path.dirname(file_path), exist_ok=True)  # Create directory if not exists
        ercc_list.save(file_path)
        ercc_list = os.path.join(app.config['UPLOAD_FOLDER'], 'ercc_list',
                                 ercc_list.filename)  # Set the file name to sirv_list
    else:
        # If no file is provided, set sirv_list to NA
        ercc_list = "NA"

    ercc_list2 = request.files.get('ercc_file2')
    print('ercc_list2', ercc_list2)
    if ercc_list2:
        # Save the file if it exists
        file_path = os.path.join(app.config['UPLOAD_FOLDER'], 'ercc_list2', ercc_list2.filename)
        os.makedirs(os.path.dirname(file_path), exist_ok=True)  # Create directory if not exists
        ercc_list2.save(file_path)
        ercc_list2 = os.path.join(app.config['UPLOAD_FOLDER'], 'ercc_list2',
                                  ercc_list2.filename)  # Set the file name to sirv_list
    else:
        # If no file is provided, set sirv_list to NA
        ercc_list2 = "NA"

    sequin_list = request.files.get('sequin_file')
    if sequin_list:
        # Save the file if it exists
        file_path = os.path.join(app.config['UPLOAD_FOLDER'], 'sequin_list', sequin_list.filename)
        os.makedirs(os.path.dirname(file_path), exist_ok=True)  # Create directory if not exists
        sequin_list.save(file_path)
        sequin_list = os.path.join(app.config['UPLOAD_FOLDER'], 'sequin_list',
                                   sequin_list.filename)  # Set the file name to sirv_list
    else:
        # If no file is provided, set sirv_list to NA
        sequin_list = "NA"

    sequin_list2 = request.files.get('sequin_file2')
    print('sequin_list2', sequin_list2)
    if sequin_list2:
        # Save the file if it exists
        file_path = os.path.join(app.config['UPLOAD_FOLDER'], 'sequin_list2', sequin_list2.filename)
        os.makedirs(os.path.dirname(file_path), exist_ok=True)  # Create directory if not exists
        sequin_list2.save(file_path)
        sequin_list2 = os.path.join(app.config['UPLOAD_FOLDER'], 'sequin_list2',
                                    sequin_list2.filename)  # Set the file name to sirv_list
    else:
        # If no file is provided, set sirv_list to NA
        sequin_list2 = "NA"

    # comparison
    selected_comparisons = request.form.getlist('comparison-dropdown')
    if not selected_comparisons:
        comparison = 'NA'
        comp_bambu = 'NA'
        comp_FLAIR = 'NA'
        comp_Lyric = 'NA'
        comp_IsoTools = 'NA'
        comp_Mandalorion = 'NA'
        comp_Iso_IB = 'NA'
        comp_FLAMES = 'NA'
        comp_IsoQuant = 'NA'
        comp_Spectra = 'NA'
        comp_TALON_LAPA = 'NA'
        comp_StringTie2 = 'NA'

    else:
        comparison = 'Custom'
        comp_bambu = 'Bambu' if 'Bambu' in selected_comparisons else 'NA'
        comp_FLAIR = 'FLAIR' if 'FLAIR' in selected_comparisons else 'NA'
        comp_Lyric = 'Lyric' if 'Lyric' in selected_comparisons else 'NA'
        comp_IsoTools = 'IsoTools' if 'IsoTools' in selected_comparisons else 'NA'
        comp_Mandalorion = 'Mandalorion' if 'Mandalorion' in selected_comparisons else 'NA'
        comp_Iso_IB = 'Iso_IB' if 'Iso_IB' in selected_comparisons else 'NA'
        comp_FLAMES = 'FLAMES' if 'FLAMES' in selected_comparisons else 'NA'
        comp_IsoQuant = 'IsoQuant' if 'IsoQuant' in selected_comparisons else 'NA'
        comp_Spectra = 'Spectra' if 'Spectra' in selected_comparisons else 'NA'
        comp_TALON_LAPA = 'TALON_LAPA' if 'TALON_LAPA' in selected_comparisons else 'NA'
        comp_StringTie2 = 'StringTie2' if 'StringTie2' in selected_comparisons else 'NA'



    print('trans_path2:', transcriptome_path2, 'anno_path2:', annotation_path2, 'reference_path2:', reference_path2, 'coverage2:', coverage2, 'coverage_dir2:', coverage_dir2)
    print('cage_path2:', cage_path2, 'quant_path2:', quant_path2, 'polyA_path2:', poly_A_path2, 'platform2:', platform2, 'lib_prep2:', library_preparation2, 'data_cat2:', data_category2)
    print('sirv_list2:', sirv_list2, 'ercc_list2:', ercc_list2, 'sequin_list2', sequin_list2)
    # Run the script in a separate thread with the additional parameters
    thread = Thread(target=run_script_process_challenge1, args=(organism, transcriptome_path, annotation_path, reference_path, coverage, coverage_dir, cage_path, quant_path, poly_A_path, platform, library_preparation, tool, data_category,
                    sirv_list, ercc_list, sequin_list, comparison, comp_bambu, comp_FLAIR, comp_Lyric, comp_IsoTools, comp_Mandalorion, comp_Iso_IB, comp_FLAMES, comp_IsoQuant, comp_Spectra, comp_TALON_LAPA, comp_StringTie2,
                                                                transcriptome_path2, annotation_path2, reference_path2, coverage2, coverage_dir2, cage_path2, quant_path2, poly_A_path2, platform2, library_preparation2, data_category2,
                                                                sirv_list2, ercc_list2, sequin_list2))

    thread.start()

    return jsonify({"status": "Challenge 1 started successfully!"})

def run_script_process_challenge1(organism, transcriptome_path, annotation_path, reference_path, coverage, coverage_dir, cage_path, quant_path, poly_A_path, platform, library_preparation, tool, data_category,
                    sirv_list, ercc_list, sequin_list, comparison, comp_bambu, comp_FLAIR, comp_Lyric, comp_IsoTools, comp_Mandalorion, comp_Iso_IB, comp_FLAMES, comp_IsoQuant, comp_Spectra, comp_TALON_LAPA, comp_StringTie2, transcriptome_path2, annotation_path2, reference_path2, coverage2, coverage_dir2,
                                  cage_path2, quant_path2, poly_A_path2, platform2, library_preparation2, data_category2, sirv_list2, ercc_list2, sequin_list2):

    print(organism, transcriptome_path, annotation_path, reference_path, coverage, coverage_dir,cage_path, quant_path,
          poly_A_path, platform, library_preparation, tool, data_category,
          sirv_list, ercc_list, sequin_list, comparison, comp_bambu, comp_FLAIR, comp_Lyric, comp_IsoTools, comp_Mandalorion, comp_Iso_IB, comp_FLAMES, comp_IsoQuant, comp_Spectra, comp_TALON_LAPA, comp_StringTie2,
          transcriptome_path2, annotation_path2, reference_path2, coverage2,
          coverage_dir2, cage_path2, quant_path2, poly_A_path2, platform2, library_preparation2, data_category2,
          sirv_list2, ercc_list2, sequin_list2)

    global progress, status_message, terminal_output

    # Define the path to the script
    script_path = "lrgasp_event2_metrics/sqanti3_lrgasp.challenge1.py"

    print('tr_path2:', transcriptome_path2)
    # Run the script
    try:

        if transcriptome_path2 == 'NA':
            print('Input received for 1 file!')
            process = subprocess.Popen(
                ['python', script_path, organism, transcriptome_path, annotation_path, reference_path, coverage, coverage_dir, cage_path, quant_path,
          poly_A_path, platform, library_preparation, tool, data_category,
          sirv_list, ercc_list, sequin_list, comparison, comp_bambu, comp_FLAIR, comp_Lyric, comp_IsoTools, comp_Mandalorion, comp_Iso_IB, comp_FLAMES, comp_IsoQuant, comp_Spectra, comp_TALON_LAPA, comp_StringTie2, transcriptome_path2,
                 annotation_path2, reference_path2, coverage2,
          coverage_dir2, cage_path2, quant_path2, poly_A_path2, platform2, library_preparation2, data_category2,
          sirv_list2, ercc_list2, sequin_list2, '--gtf'],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True
            )
        else:
            print('Input received for 2 files!')
            process = subprocess.Popen(
                ['python', script_path, organism, transcriptome_path, annotation_path, reference_path, coverage,
                 coverage_dir, cage_path, quant_path,
                 poly_A_path, platform, library_preparation, tool, data_category,
                 sirv_list, ercc_list, sequin_list, comparison, comp_bambu, comp_FLAIR, comp_Lyric, comp_IsoTools, comp_Mandalorion, comp_Iso_IB, comp_FLAMES, comp_IsoQuant, comp_Spectra, comp_TALON_LAPA, comp_StringTie2, transcriptome_path2,
                 annotation_path2, reference_path2, coverage2,
                 coverage_dir2, cage_path2, quant_path2, poly_A_path2, platform2, library_preparation2, data_category2,
                 sirv_list2, ercc_list2, sequin_list2, '--dataset2', '--gtf'],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True
            )

        for line in process.stdout:
            line = line.strip()
            print(line)
            terminal_output.append(line)

            if 'PROGRESS:' in line:
                progress = int(line.strip().split('PROGRESS:')[1])
                print(f"Progress updated to {progress}%")
            elif 'STATUS:' in line:
                status_message = line.strip().split('STATUS:')[1].strip()
                print(f"Status message updated to: {status_message}")

        process.wait()
        progress = 100
        status_message = 'Script Finished'

        # Move the generated report
        source = glob.glob("sqanti_results/results_file1/*.html")[0]
        destination = "static/report.html"
        shutil.move(source, destination)

    except Exception as e:
        print(f"Error running challenge 1: {e}")
        status_message = "Error running script"
        progress = 0  # Ensure progress reaches 100 on error

@app.route("/challenge3", methods=['GET', 'POST'])
def challenge3():
    return render_template("challenge_3_new.html")

@app.route("/run_script_challenge3", methods=["POST"])
def run_script_challenge3():
    global progress, status_message

    # Reset progress and status message
    progress = 0
    status_message = "Starting Script..."

    def replace_directory(directory_path):
        if os.path.exists(directory_path):
            shutil.rmtree(directory_path)
        os.makedirs(directory_path)

    replace_directory("uploads")
    replace_directory("uploads/transcriptome_file1")
    replace_directory("uploads/transcriptome_file2")
    replace_directory("sqanti_results/results_file1")
    replace_directory("sqanti_results/results_file2")
    replace_directory("uploads/coverage_files")
    replace_directory("uploads/coverage_files2")

    # Get file and form data from the request
    file = request.files.get('file')
    transcriptome_path = os.path.join(app.config['UPLOAD_FOLDER'], 'transcriptome_file1', file.filename)
    file.save(transcriptome_path)

    # meta data
    organism = request.form.get('organism')
    platform = request.form.get('platform')
    tool = request.form.get('tool')
    library_preparation = request.form.get('library_preparation')
    data_category = request.form.get('data_category')

    # comparison
    selected_comparisons = request.form.getlist('comparison-dropdown')

    if not selected_comparisons:
        comparison = 'NA'
        comp_bambu = 'NA'
        comp_RNA_Bloom = 'NA'
        comp_rnaSPAdes = 'NA'
        comp_StringTie2_IsoQuant = 'NA'
    else:
        comparison = 'Custom'
        comp_bambu = 'Bambu' if 'Bambu' in selected_comparisons else 'NA'
        comp_RNA_Bloom = 'RNA_Bloom' if 'RNA_Bloom' in selected_comparisons else 'NA'
        comp_rnaSPAdes = 'rnaSPAdes' if 'rnaSPAdes' in selected_comparisons else 'NA'
        comp_StringTie2_IsoQuant = 'StringTie2_IsoQuant' if 'StringTie2_IsoQuant' in selected_comparisons else 'NA'


    annotation_file = request.form.get('annotation')
    if annotation_file == 'custom':
        annotation_file = request.files.get('annotation_file')
        annotation_path = os.path.join(app.config['UPLOAD_FOLDER'], 'transcriptome_file1', annotation_file.filename)
        annotation_file.save(annotation_path)
    elif annotation_file == 'default':
        annotation_path = 'LRGASP_DATA'
    else:
        annotation_path = 'NA'

    reference_file = request.files.get('reference_file')
    if reference_file and reference_file.filename != '':
        reference_path = os.path.join(app.config['UPLOAD_FOLDER'], 'transcriptome_file1', reference_file.filename)
        reference_file.save(reference_path)
    else:
        reference_path = 'LRGASP_DATA'

    coverage = request.form.get('coverage_directory')
    if coverage == 'custom':
        coverage_files = request.files.getlist('coverage_dir')
        coverage_dir = coverage_files[0].filename.split('/')[0]
        replace_directory(os.path.join("uploads/coverage_files", coverage_dir))
        for cov_file in coverage_files:
            coverage_path = os.path.join(app.config['UPLOAD_FOLDER'], 'coverage_files', cov_file.filename)
            cov_file.save(coverage_path)
    elif coverage == 'default':
        coverage = 'LRGASP_DATA'
        coverage_dir = 'LRGASP_DATA'
    else:
        coverage = 'NA'
        coverage_dir = 'NA'

    # All variables for the 2nd dataset
    transcriptome_file2 = request.files.get('file-2')
    if transcriptome_file2.filename != '':
        transcriptome_path2 = os.path.join(app.config['UPLOAD_FOLDER'], 'transcriptome_file2', transcriptome_file2.filename)
        transcriptome_file2.save(transcriptome_path2)
    else:
        transcriptome_path2 = 'NA'

    platform2 = request.form.get('platform2')
    library_preparation2 = request.form.get('library_preparation2')
    data_category2 = request.form.get('data_category2')

    annotation_file2 = request.form.get('annotation-2')
    if annotation_file2 == 'custom':
        annotation_file2 = request.files.get('annotation_file2')
        annotation_path2 = os.path.join(app.config['UPLOAD_FOLDER'], 'transcriptome_file2', annotation_file2.filename)
        annotation_file2.save(annotation_path2)
    elif annotation_file2 == 'default':
        annotation_path2 = 'LRGASP_DATA'
    else:
        annotation_path2 = 'NA'

    reference_file_2 = request.form.get('reference-2')
    if reference_file_2 == 'custom':
        reference_file_2 = request.files.get('reference_file2')
        reference_path2 = os.path.join(app.config['UPLOAD_FOLDER'], 'transcriptome_file2', reference_file_2.filename)
        reference_file_2.save(reference_path2)
    else:
        reference_path2 = 'LRGASP_DATA'

    coverage2 = request.form.get('coverage_directory2')
    if coverage2== 'custom':
        coverage_files2 = request.files.getlist('coverage_dir2')
        coverage_dir2 = coverage_files2[0].filename.split('/')[0]
        replace_directory(os.path.join("uploads/coverage_files2", coverage_dir2))
        for cov_file in coverage_files2:
            coverage_path = os.path.join(app.config['UPLOAD_FOLDER'], 'coverage_files2', cov_file.filename)
            cov_file.save(coverage_path)
    elif coverage2 == 'default':
        coverage2 = 'LRGASP_DATA'
        coverage_dir2 = 'LRGASP_DATA'
    else:
        coverage2 = 'NA'
        coverage_dir2 = 'NA'

    sirv_list = request.files.get('sirv_file1')
    if sirv_list:
        # Save the file if it exists
        file_path = os.path.join(app.config['UPLOAD_FOLDER'], 'sirv_list', sirv_list.filename)
        os.makedirs(os.path.dirname(file_path), exist_ok=True)  # Create directory if not exists
        sirv_list.save(file_path)
        sirv_list = os.path.join(app.config['UPLOAD_FOLDER'], 'sirv_list', sirv_list.filename)
    else:
        # If no file is provided, set sirv_list to NA
        sirv_list = "NA"

    ercc_list = request.files.get('ercc_file1')
    if ercc_list:
        # Save the file if it exists
        file_path = os.path.join(app.config['UPLOAD_FOLDER'], 'ercc_list', ercc_list.filename)
        os.makedirs(os.path.dirname(file_path), exist_ok=True)  # Create directory if not exists
        ercc_list.save(file_path)
        ercc_list = os.path.join(app.config['UPLOAD_FOLDER'], 'ercc_list', ercc_list.filename)  # Set the file name to sirv_list
    else:
        # If no file is provided, set sirv_list to NA
        ercc_list = "NA"

    sequin_list = request.files.get('sequin_file1')
    if sequin_list:
        # Save the file if it exists
        file_path = os.path.join(app.config['UPLOAD_FOLDER'], 'sequin_list', sequin_list.filename)
        os.makedirs(os.path.dirname(file_path), exist_ok=True)  # Create directory if not exists
        sequin_list.save(file_path)
        sequin_list = os.path.join(app.config['UPLOAD_FOLDER'], 'sequin_list', sequin_list.filename)  # Set the file name to sirv_list
    else:
        # If no file is provided, set sirv_list to NA
        sequin_list = "NA"

    sirv_list2 = request.files.get('sirv_file2')
    if sirv_list2:
        # Save the file if it exists
        file_path = os.path.join(app.config['UPLOAD_FOLDER'], 'sirv_list2', sirv_list2.filename)
        os.makedirs(os.path.dirname(file_path), exist_ok=True)  # Create directory if not exists
        sirv_list2.save(file_path)
        sirv_list2 = os.path.join(app.config['UPLOAD_FOLDER'], 'sirv_list2', sirv_list2.filename)
    else:
        # If no file is provided, set sirv_list to NA
        sirv_list2 = "NA"

    ercc_list2 = request.files.get('ercc_file2')
    if ercc_list2:
        # Save the file if it exists
        file_path = os.path.join(app.config['UPLOAD_FOLDER'], 'ercc_list2', ercc_list2.filename)
        os.makedirs(os.path.dirname(file_path), exist_ok=True)  # Create directory if not exists
        ercc_list2.save(file_path)
        ercc_list2 = os.path.join(app.config['UPLOAD_FOLDER'], 'ercc_list2',
                                 ercc_list2.filename)  # Set the file name to sirv_list
    else:
        # If no file is provided, set sirv_list to NA
        ercc_list2 = "NA"

    sequin_list2 = request.files.get('sequin_file2')
    if sequin_list2:
        # Save the file if it exists
        file_path = os.path.join(app.config['UPLOAD_FOLDER'], 'sequin_list2', sequin_list2.filename)
        os.makedirs(os.path.dirname(file_path), exist_ok=True)  # Create directory if not exists
        sequin_list2.save(file_path)
        sequin_list2 = os.path.join(app.config['UPLOAD_FOLDER'], 'sequin_list2',
                                   sequin_list2.filename)  # Set the file name to sirv_list
    else:
        # If no file is provided, set sirv_list to NA
        sequin_list2 = "NA"

    # 2nd dataset provided?
    dataset2 = request.form.get('dataset2')
    if dataset2 != 'none':
        dataset2 = True
    else:
        dataset2 = False

    # Run the script in a separate thread with the additional parameters
    thread = Thread(target=run_script_process_challenge3, args=(transcriptome_path, organism, platform, library_preparation, tool, data_category, annotation_path, reference_path, coverage, coverage_dir,
                                                     transcriptome_path2, platform2, library_preparation2, data_category2, annotation_path2, reference_path2, coverage2, coverage_dir2,
                                                     comparison, comp_bambu, comp_RNA_Bloom, comp_rnaSPAdes, comp_StringTie2_IsoQuant,
                                                     sirv_list, ercc_list, sequin_list, sirv_list2, ercc_list2, sequin_list2, dataset2))
    thread.start()

    return jsonify({"status": "Challenge 3 started successfully!"})


def run_script_process_challenge3(file_path, organism, platform, library_preparation, tool, data_category, annotation_path, reference_path, coverage, coverage_dir,
                       file_path_2, platform2, library_preparation2, data_category2, annotation_path_2, reference_path_2, coverage2, coverage_dir2,
                       comparison, comp_bambu, comp_RNA_Bloom, comp_rnaSPAdes, comp_StringTie2_IsoQuant,
                       sirv_list, ercc_list, sequin_list, sirv_list2, ercc_list2, sequin_list2, dataset2):

    global progress, status_message, terminal_output

    terminal_output = []

    # Define the path to the script
    script_path = "lrgasp_event2_metrics/sqanti3_lrgasp.challenge3.py"

    # Run the script
    try:
        if dataset2 == False:
            print('########### 1 DATASET ############')
            process = subprocess.Popen(
                ['python', script_path, file_path, organism, platform, library_preparation, tool, data_category, annotation_path, reference_path, coverage, coverage_dir,
                                file_path_2, platform2, library_preparation2, data_category2, annotation_path_2, reference_path_2, coverage2, coverage_dir2,
                                comparison, comp_bambu, comp_RNA_Bloom, comp_rnaSPAdes, comp_StringTie2_IsoQuant,
                                sirv_list, ercc_list, sequin_list, sirv_list2, ercc_list2, sequin_list2],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True
            )
        elif dataset2 == True:
            print('########### 2 DATASETS ############')
            process = subprocess.Popen(
                ['python', script_path, file_path, organism, platform, library_preparation, tool, data_category, annotation_path, reference_path, coverage, coverage_dir,
                                file_path_2, platform2, library_preparation2, data_category2, annotation_path_2, reference_path_2, coverage2, coverage_dir2,
                                comparison, comp_bambu, comp_RNA_Bloom, comp_rnaSPAdes, comp_StringTie2_IsoQuant,
                                sirv_list, ercc_list, sequin_list, sirv_list2, ercc_list2, sequin_list2, '--dataset2'],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True
            )
        else:
            print('ERROR: Could not detect if 1 or 2 datasets are provided! Exiting...')


        for line in process.stdout:
            line = line.strip()
            print(line)
            terminal_output.append(line)

            if 'PROGRESS:' in line:
                progress = int(line.strip().split('PROGRESS:')[1])
                print(f"Progress updated to {progress}%")
            elif 'STATUS:' in line:
                status_message = line.strip().split('STATUS:')[1].strip()
                print(f"Status message updated to: {status_message}")

        process.wait()
        progress = 100
        status_message = 'Script Finished'

        # Move the generated report
        source = glob.glob("sqanti_results/results_file1/*.html")[0]
        destination = "static/report.html"
        shutil.move(source, destination)

    except Exception as e:
        print(f"Error running script: {e}")
        status_message = "Error running script"
        progress = 0


@app.route("/terminal_stream")
def terminal_stream():

    def generate():
        last_index = 0
        while True:
            # If new lines have been appended, yield them
            if last_index < len(terminal_output):
                new_lines = terminal_output[last_index:]
                for line in new_lines:
                    yield f"data: {line}\n\n"
                last_index = len(terminal_output)
            time.sleep(0.5)
    return Response(generate(), mimetype="text/event-stream")

@app.route("/progress_stream")
def progress_stream():
    """Stream progress updates to the client."""
    def generate():
        while True:
            global progress
            global status_message
            time.sleep(0.5)  # Check progress every 500ms
            data = json.dumps({'progress': progress, 'message': status_message})
            yield f"data: {data}\n\n"
            if progress >= 100:
                break
    return Response(generate(), mimetype="text/event-stream")



@app.route('/challenge_results')
def challenge3_results():
    return render_template("show_report.html")

@app.route('/download-report')
def download_report():
    return send_from_directory(os.path.join(app.root_path, 'static'), 'report.html', as_attachment=True)


if __name__ == '__main__':
    app.run(debug=True)
