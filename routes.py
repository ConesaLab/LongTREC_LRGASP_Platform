from flask import render_template, redirect, url_for, flash, request, jsonify, send_file, current_app
from flask_login import login_user, logout_user, login_required, current_user
from app import app, db
from models import User, Submission, Job
from forms import RegistrationForm, LoginForm, SubmissionForm, Challenge3Form, Challenge1Form
from werkzeug.security import generate_password_hash, check_password_hash
from werkzeug.utils import secure_filename
import os
import sys
import subprocess
import glob
import shutil
from threading import Thread
from datetime import datetime

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
        hashed_password = generate_password_hash(form.password.data)
        user = User(username=form.username.data, email=form.email.data, password=hashed_password)
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
        if user and check_password_hash(user.password, form.password.data):
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

@app.route('/challenge3', methods=['GET', 'POST'])
def challenge3():
    form = Challenge3Form()
    return render_template('challenge_3_new.html', title='Challenge 3', form=form)

@app.route('/challenge1', methods=['GET', 'POST'])
def challenge1():
    form = Challenge1Form()
    return render_template('challenge1_final.html', title='Challenge 1', form=form)

@app.route('/run_script_challenge3', methods=['POST'])
def run_script_challenge3():
    try:
        email = request.form.get('email')
        if not email:
            return jsonify({'error': 'Email is required'}), 400
        
        job = Job(
            user_id=current_user.id if current_user.is_authenticated else None,
            email=email,
            challenge_type='challenge3',
            status='pending'
        )
        db.session.add(job)
        db.session.commit()
        
        job_id = job.id
        
        from utils import send_start_email
        send_start_email(job_id)
        
        # Helper function to replace directories
        def replace_directory(directory_path):
            if os.path.exists(directory_path):
                shutil.rmtree(directory_path)
            os.makedirs(directory_path)

        job_results_dir = f"sqanti_results/{job_id}"
        results_file1_dir = f"{job_results_dir}/results_file1"
        results_file2_dir = f"{job_results_dir}/results_file2"
        
        replace_directory(results_file1_dir)
        replace_directory(results_file2_dir)
        replace_directory(f"{results_file1_dir}/coverage_files")
        replace_directory(f"{results_file2_dir}/coverage_files")

        # Process transcriptome file 1
        file = request.files.get('file', default='')
        if file.filename == '':
            return jsonify({'error': 'No valid transcriptome file uploaded'}), 400
        transcriptome_path = os.path.join(results_file1_dir, file.filename)
        file.save(transcriptome_path)

        # Get metadata
        organism = request.form.get('organism')
        platform = request.form.get('platform')
        tool = request.form.get('tool')
        library_preparation = request.form.get('library_preparation')
        data_category = request.form.get('data_category')

        # Process comparison options
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

        # Process annotation file
        annotation_file = request.form.get('annotation')
        if annotation_file == 'custom':
            annotation_file = request.files.get('annotation_file')
            annotation_path = os.path.join(results_file1_dir, annotation_file.filename)
            annotation_file.save(annotation_path)
        elif annotation_file == 'default':
            annotation_path = 'LRGASP_DATA'
        else:
            annotation_path = 'NA'

        # Process reference file
        reference_file = request.files.get('reference_file')
        if reference_file and reference_file.filename != '':
            reference_path = os.path.join(results_file1_dir, reference_file.filename)
            reference_file.save(reference_path)
        else:
            reference_path = 'LRGASP_DATA'

        # Process coverage directory
        coverage = request.form.get('coverage_directory')
        if coverage == 'custom':
            coverage_files = request.files.getlist('coverage_dir')
            coverage_dir = coverage_files[0].filename.split('/')[0]
            replace_directory(os.path.join(results_file1_dir, "coverage_files", coverage_dir))
            for cov_file in coverage_files:
                coverage_path = os.path.join(results_file1_dir, 'coverage_files', cov_file.filename)
                cov_file.save(coverage_path)
            coverage_dir = os.path.join(results_file1_dir, "coverage_files", coverage_dir)
        elif coverage == 'default':
            coverage = 'LRGASP_DATA'
            coverage_dir = 'LRGASP_DATA'
        else:
            coverage = 'NA'
            coverage_dir = 'NA'

        # Check if dataset2 is provided
        dataset2 = request.form.get('dataset2')
        if dataset2 != 'none':
            dataset2 = True
        else:
            dataset2 = False

        # Process transcriptome file 2
        transcriptome_file2 = request.files.get('file-2')
        if transcriptome_file2 and transcriptome_file2.filename != '':
            file_path_2 = os.path.join(results_file2_dir, transcriptome_file2.filename)
            transcriptome_file2.save(file_path_2)
        else:
            file_path_2 = 'NA'

        # Get metadata for dataset 2
        platform2 = request.form.get('platform2', 'NA')
        library_preparation2 = request.form.get('library_preparation2', 'NA')
        data_category2 = request.form.get('data_category2', 'NA')

        # Process annotation file 2
        annotation_file_2 = request.form.get('annotation-2')
        if annotation_file_2 == 'custom':
            annotation_file_2 = request.files.get('annotation_file2')
            annotation_path_2 = os.path.join(results_file2_dir, annotation_file_2.filename)
            if dataset2:
                annotation_file_2.save(annotation_path_2)
        elif annotation_file_2 == 'default':
            annotation_path_2 = 'LRGASP_DATA'
        else:
            annotation_path_2 = 'NA'

        # Process reference file 2
        reference_file_2 = request.files.get('reference-2')
        if reference_file_2 and reference_file_2.filename != '':
            reference_path_2 = os.path.join(results_file2_dir, reference_file_2.filename)
            if dataset2:
                reference_file_2.save(reference_path_2)
        else:
            reference_path_2 = 'LRGASP_DATA'

        # Process coverage directory 2
        coverage2 = request.form.get('coverage_directory2')
        if coverage2 == 'custom':
            coverage_files2 = request.files.getlist('coverage_dir2')
            coverage_dir2 = coverage_files2[0].filename.split('/')[0]
            replace_directory(os.path.join(results_file2_dir, "coverage_files", coverage_dir2))
            for cov_file in coverage_files2:
                coverage_path = os.path.join(results_file2_dir, 'coverage_files', cov_file.filename)
                if dataset2:
                    cov_file.save(coverage_path)
            coverage_dir2 = os.path.join(results_file2_dir, "coverage_files", coverage_dir2)
        elif coverage2 == 'default':
            coverage2 = 'LRGASP_DATA'
            coverage_dir2 = 'LRGASP_DATA'
        else:
            coverage2 = 'NA'
            coverage_dir2 = 'NA'

        # Process SIRV list files
        sirv_list = request.files.get('sirv_file1')
        if sirv_list:
            file_path = os.path.join(results_file1_dir, sirv_list.filename)
            sirv_list.save(file_path)
            sirv_list = file_path
        else:
            sirv_list = "NA"

        sirv_list2 = request.files.get('sirv_file2')
        if sirv_list2:
            file_path = os.path.join(results_file2_dir, sirv_list2.filename)
            sirv_list2.save(file_path)
            sirv_list2 = file_path
        else:
            sirv_list2 = "NA"

        # Process ERCC list files
        ercc_list = request.files.get('ercc_file1')
        if ercc_list:
            file_path = os.path.join(results_file1_dir, ercc_list.filename)
            ercc_list.save(file_path)
            ercc_list = file_path
        else:
            ercc_list = "NA"

        ercc_list2 = request.files.get('ercc_file2')
        if ercc_list2:
            file_path = os.path.join(results_file2_dir, ercc_list2.filename)
            ercc_list2.save(file_path)
            ercc_list2 = file_path
        else:
            ercc_list2 = "NA"

        # Process Sequin list files
        sequin_list = request.files.get('sequin_file1')
        if sequin_list:
            file_path = os.path.join(results_file1_dir, sequin_list.filename)
            sequin_list.save(file_path)
            sequin_list = file_path
        else:
            sequin_list = "NA"

        sequin_list2 = request.files.get('sequin_file2')
        if sequin_list2:
            file_path = os.path.join(results_file2_dir, sequin_list2.filename)
            sequin_list2.save(file_path)
            sequin_list2 = file_path
        else:
            sequin_list2 = "NA"
        
        thread = Thread(target=run_script_process_challenge3, args=(
            job.id, job_results_dir, transcriptome_path, organism, platform, library_preparation, tool, data_category,
            annotation_path, reference_path, coverage, coverage_dir,
            file_path_2, platform2, library_preparation2, data_category2,
            annotation_path_2, reference_path_2, coverage2, coverage_dir2,
            comparison, comp_bambu, comp_RNA_Bloom, comp_rnaSPAdes, comp_StringTie2_IsoQuant,
            sirv_list, ercc_list, sequin_list, sirv_list2, ercc_list2, sequin_list2, dataset2
        ))
        thread.daemon = True
        thread.start()
        
        return jsonify({
            'success': True,
            'message': 'Your analysis has been started. You will receive an email when it completes.',
            'job_id': job.id
        })
        
    except Exception as e:
        print(f"Error starting job: {str(e)}", file=sys.stderr)
        return jsonify({'error': str(e)}), 500

def run_script_process_challenge3(job_id, job_results_dir, transcriptome_path, organism, platform, library_preparation, tool, data_category,
                                  annotation_path, reference_path, coverage, coverage_dir,
                                  file_path_2, platform2, library_preparation2, data_category2,
                                  annotation_path_2, reference_path_2, coverage2, coverage_dir2,
                                  comparison, comp_bambu, comp_RNA_Bloom, comp_rnaSPAdes, comp_StringTie2_IsoQuant,
                                  sirv_list, ercc_list, sequin_list, sirv_list2, ercc_list2, sequin_list2, dataset2):
    """Background process for Challenge 3"""
    with app.app_context():
        job = Job.query.get(job_id)
        if not job:
            return
        
        logs_dir = "logs"
        os.makedirs(logs_dir, exist_ok=True)
        log_file_path = os.path.join(logs_dir, f"{job_id}.log")
        
        try:
            with open(log_file_path, 'w') as log_file:
                log_file.write(f"=== Challenge 3 Analysis Started ===\n")
                log_file.write(f"Job ID: {job_id}\n")
                log_file.write(f"Results Directory: {job_results_dir}\n")
                log_file.write(f"Started at: {datetime.utcnow().isoformat()}\n")
                log_file.write(f"Organism: {organism}\n")
                log_file.write(f"Platform: {platform}\n")
                log_file.write(f"Tool: {tool}\n")
                log_file.write(f"Library Preparation: {library_preparation}\n")
                log_file.write(f"Data Category: {data_category}\n")
                log_file.write(f"Dataset 2: {dataset2}\n")
                log_file.write(f"{'='*50}\n\n")
                log_file.flush()
                
                job.status = 'running'
                job.started_at = datetime.utcnow()
                db.session.commit()
                
                script_path = "lrgasp_event2_metrics/sqanti3_lrgasp.challenge3.py"
                
                # Build command arguments
                cmd_args = [
                    'python', script_path,  str(job_id),
                    transcriptome_path, organism, platform, library_preparation, tool, data_category,
                    annotation_path, reference_path, coverage, coverage_dir,
                    file_path_2, platform2, library_preparation2, data_category2,
                    annotation_path_2, reference_path_2, coverage2, coverage_dir2,
                    comparison, comp_bambu, comp_RNA_Bloom, comp_rnaSPAdes, comp_StringTie2_IsoQuant,
                    sirv_list, ercc_list, sequin_list, sirv_list2, ercc_list2, sequin_list2
                ]
                
                # Add --dataset2 flag if second dataset is provided
                if dataset2:
                    cmd_args.append('--dataset2')
                
                log_file.write(f"Executing command:\n{' '.join(cmd_args)}\n\n")
                log_file.write(f"{'='*50}\n")
                log_file.write(f"Script Output:\n")
                log_file.write(f"{'='*50}\n\n")
                log_file.flush()
                
                process = subprocess.Popen(
                    cmd_args,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    cwd=os.getcwd()
                )
                
                for line in process.stdout:
                    line = line.rstrip()
                    print(line)
                    log_file.write(line + '\n')
                    log_file.flush()
                
                process.wait()
                
                log_file.write(f"\n{'='*50}\n")
                log_file.write(f"Script finished with return code: {process.returncode}\n")
                log_file.flush()
                
                if process.returncode != 0:
                    raise Exception(f"Script exited with code {process.returncode}")
                
                report_pattern = f"{job_results_dir}/results_file1/*.html"
                reports = glob.glob(report_pattern)
                if not reports:
                    raise Exception(f"No HTML report found in {job_results_dir}/results_file1/")
                
                source = reports[0]
                destination = f"static/report_{job_id}.html"
                shutil.copy(source, destination)
                
                log_file.write(f"Report copied from {source} to: {destination}\n")
                log_file.flush()
                
                job.status = 'completed'
                job.completed_at = datetime.utcnow()
                job.result_path = destination
                db.session.commit()
                
                log_file.write(f"\n{'='*50}\n")
                log_file.write(f"=== Challenge 3 Analysis Completed Successfully ===\n")
                log_file.write(f"Completed at: {datetime.utcnow().isoformat()}\n")
                log_file.write(f"{'='*50}\n")
                
                from utils import send_completion_email
                send_completion_email(job.id)
            
        except Exception as e:
            error_msg = f"Error in job {job_id}: {str(e)}"
            print(error_msg, file=sys.stderr)
            
            with open(log_file_path, 'a') as log_file:
                log_file.write(f"\n{'='*50}\n")
                log_file.write(f"=== ERROR ===\n")
                log_file.write(f"{error_msg}\n")
                log_file.write(f"{'='*50}\n")
            
            job.status = 'failed'
            job.completed_at = datetime.utcnow()
            job.error_message = str(e)
            db.session.commit()
            
            from utils import send_failure_email
            send_failure_email(job.id)

@app.route('/run_script_challenge1', methods=['POST'])
def run_script_challenge1():
    try:
        email = request.form.get('email')
        if not email:
            return jsonify({'error': 'Email is required'}), 400
        
        job = Job(
            user_id=current_user.id if current_user.is_authenticated else None,
            email=email,
            challenge_type='challenge1',
            status='pending'
        )
        db.session.add(job)
        db.session.commit()
        
        job_id = job.id
        
        from utils import send_start_email
        send_start_email(job_id)
        
        # Helper function to replace directories
        def replace_directory(directory_path):
            if os.path.exists(directory_path):
                shutil.rmtree(directory_path)
            os.makedirs(directory_path)

        job_upload_dir = f"uploads/job_{job_id}"
        job_results_dir = f"sqanti_results/{job_id}"
        
        replace_directory(job_upload_dir)
        replace_directory(f"{job_upload_dir}/transcriptome_file1")
        replace_directory(f"{job_upload_dir}/transcriptome_file2")
        replace_directory(f"{job_results_dir}/results_file1")
        replace_directory(f"{job_results_dir}/results_file2")
        replace_directory(f"{job_upload_dir}/coverage_files")
        replace_directory(f"{job_upload_dir}/coverage_files2")

        organism = request.form.get('organism')

        # Process transcriptome file 1
        transcriptome_file = request.files.get('file')
        transcriptome_path = os.path.join(job_upload_dir, 'transcriptome_file1', transcriptome_file.filename)
        transcriptome_file.save(transcriptome_path)

        # Process transcriptome file 2
        transcriptome_file2 = request.files.get('file-2')
        if transcriptome_file2 and transcriptome_file2.filename != '':
            transcriptome_path2 = os.path.join(job_upload_dir, 'transcriptome_file2', transcriptome_file2.filename)
            transcriptome_file2.save(transcriptome_path2)
        else:
            transcriptome_path2 = 'NA'


        # Process annotation file 1
        annotation_file = request.form.get('annotation')
        if annotation_file == 'default':
            annotation_path = 'LRGASP_DATA'
        else:
            annotation_file = request.files.get('annotation_file')
            annotation_path = os.path.join(job_upload_dir, 'transcriptome_file1', annotation_file.filename)
            annotation_file.save(annotation_path)

        # Process annotation file 2
        annotation_file2 = request.form.get('annotation2')
        if annotation_file2 != None:
            if annotation_file2 == 'default':
                annotation_path2 = 'LRGASP_DATA'
            else:
                annotation_file2 = request.files.get('annotation_file-2')
                annotation_path2 = os.path.join(job_upload_dir, 'transcriptome_file2', annotation_file2.filename)
                annotation_file2.save(annotation_path2)
        else:
            annotation_path2 = 'NA'


        # Process reference file 1
        reference_file = request.files.get('reference_file')
        if reference_file and reference_file.filename != '':
            reference_path = os.path.join(job_upload_dir, 'transcriptome_file1', reference_file.filename)
            reference_file.save(reference_path)
        else:
            reference_path = 'LRGASP_DATA'

        # Process reference file 2
        reference_file2 = request.files.get('reference_file-2')
        if reference_file2 != None:
            if reference_file2 and reference_file2.filename != '':
                reference_path2 = os.path.join(job_upload_dir, 'transcriptome_file2', reference_file2.filename)
                reference_file2.save(reference_path2)
            else:
                reference_path2 = 'LRGASP_DATA'
        else:
            reference_path2 = 'NA'


        # Process coverage directory 1
        coverage = request.form.get('coverage_directory')
        if coverage == 'custom':
            coverage_files = request.files.getlist('coverage_directory[]')
            coverage_dir = coverage_files[0].filename.split('/')[0]
            replace_directory(os.path.join(job_upload_dir, "coverage_files", coverage_dir))
            for cov_file in coverage_files:
                coverage_path = os.path.join(job_upload_dir, 'coverage_files', cov_file.filename)
                cov_file.save(coverage_path)
        else:
            coverage_dir = 'NA'

        # Process coverage directory 2
        coverage2 = request.form.get('coverage_directory-2')
        if coverage2 == 'custom':
            coverage_files2 = request.files.getlist('coverage_directory_2[]')
            coverage_dir2 = coverage_files2[0].filename.split('/')[0]
            replace_directory(os.path.join(job_upload_dir, "coverage_files2", coverage_dir2))
            for cov_file in coverage_files2:
                coverage_path2 = os.path.join(job_upload_dir, 'coverage_files2', cov_file.filename)
                cov_file.save(coverage_path2)
        else:
            coverage_dir2 = 'NA'
            coverage2 = 'NA'


        # Process CAGE file 1
        cage_file = request.form.get('cage')
        if cage_file == 'default':
            cage_path = 'LRGASP_DATA'
        elif cage_file == 'reference':
            cage_path = 'reference'
        else:
            cage_file = request.files.get('cage_file')
            cage_path = os.path.join(job_upload_dir, 'transcriptome_file1', cage_file.filename)
            cage_file.save(cage_path)

        # Process CAGE file 2
        cage_file2 = request.form.get('cage2')
        if cage_file2 != None:
            if cage_file2 == 'default':
                cage_path2 = 'LRGASP_DATA'
            elif cage_file2 == 'reference':
                cage_path2 = 'reference'
            else:
                cage_file2 = request.files.get('cage_file2')
                cage_path2 = os.path.join(job_upload_dir, 'transcriptome_file2', cage_file2.filename)
                cage_file2.save(cage_path2)
        else:
            cage_path2 = 'NA'


        # Process quant file 1
        quant_file = request.form.get('quant')
        if quant_file == 'default':
            quant_path = 'LRGASP_DATA'
        elif quant_file == 'reference':
            quant_path = 'reference'
        else:
            quant_file = request.files.get('quant_file')
            quant_path = os.path.join(job_upload_dir, 'transcriptome_file1', quant_file.filename)
            quant_file.save(quant_path)

        # Process quant file 2
        quant_file2 = request.form.get('quant2')
        if quant_file2 != None:
            if quant_file2 == 'default':
                quant_path2 = 'LRGASP_DATA'
            elif quant_file2 == 'reference':
                quant_path2 = 'reference'
            else:
                quant_file2 = request.files.get('quant_file2')
                quant_path2 = os.path.join(job_upload_dir, 'transcriptome_file2', quant_file2.filename)
                quant_file2.save(quant_path2)
        else:
            quant_path2 = 'NA'


        # Process polyA file 1
        poly_A_file = request.form.get('polyA')
        if poly_A_file == 'default':
            poly_A_path = 'LRGASP_DATA'
        elif poly_A_file == 'NA':
            poly_A_path = 'NA'
        else:
            poly_A_file = request.files.get('polyA_file')
            poly_A_path = os.path.join(job_upload_dir, 'transcriptome_file1', poly_A_file.filename)
            poly_A_file.save(poly_A_path)

        # Process polyA file 2
        poly_A_file2 = request.form.get('polyA2')
        if poly_A_file2 != None:
            if poly_A_file2 == 'default':
                poly_A_path2 = 'LRGASP_DATA'
            elif poly_A_file2 == 'NA':
                poly_A_path2 = 'NA'
            else:
                poly_A_file2 = request.files.get('polyA_file2')
                poly_A_path2 = os.path.join(job_upload_dir, 'transcriptome_file2', poly_A_file2.filename)
                poly_A_file2.save(poly_A_path2)
        else:
            poly_A_path2 = 'NA'


        # Get metadata
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


        # Process SIRV list files
        sirv_list = request.files.get('sirv_file')
        if sirv_list:
            file_path = os.path.join(job_upload_dir, 'sirv_list', sirv_list.filename)
            os.makedirs(os.path.dirname(file_path), exist_ok=True)
            sirv_list.save(file_path)
            sirv_list = file_path
        else:
            sirv_list = "NA"

        sirv_list2 = request.files.get('sirv_file2')
        if sirv_list2:
            file_path = os.path.join(job_upload_dir, 'sirv_list2', sirv_list2.filename)
            os.makedirs(os.path.dirname(file_path), exist_ok=True)
            sirv_list2.save(file_path)
            sirv_list2 = file_path
        else:
            sirv_list2 = "NA"


        # Process ERCC list files
        ercc_list = request.files.get('ercc_file')
        if ercc_list:
            file_path = os.path.join(job_upload_dir, 'ercc_list', ercc_list.filename)
            os.makedirs(os.path.dirname(file_path), exist_ok=True)
            ercc_list.save(file_path)
            ercc_list = file_path
        else:
            ercc_list = "NA"

        ercc_list2 = request.files.get('ercc_file2')
        if ercc_list2:
            file_path = os.path.join(job_upload_dir, 'ercc_list2', ercc_list2.filename)
            os.makedirs(os.path.dirname(file_path), exist_ok=True)
            ercc_list2.save(file_path)
            ercc_list2 = file_path
        else:
            ercc_list2 = "NA"


        # Process Sequin list files
        sequin_list = request.files.get('sequin_file')
        if sequin_list:
            file_path = os.path.join(job_upload_dir, 'sequin_list', sequin_list.filename)
            os.makedirs(os.path.dirname(file_path), exist_ok=True)
            sequin_list.save(file_path)
            sequin_list = file_path
        else:
            sequin_list = "NA"

        sequin_list2 = request.files.get('sequin_file2')
        if sequin_list2:
            file_path = os.path.join(job_upload_dir, 'sequin_list2', sequin_list2.filename)
            os.makedirs(os.path.dirname(file_path), exist_ok=True)
            sequin_list2.save(file_path)
            sequin_list2 = file_path
        else:
            sequin_list2 = "NA"


        # Process comparison options
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
        
        thread = Thread(target=run_script_process_challenge1, args=(
            job.id, job_results_dir, organism, transcriptome_path, annotation_path, reference_path,
            coverage, coverage_dir, cage_path, quant_path, poly_A_path, platform,
            library_preparation, tool, data_category, sirv_list, ercc_list, sequin_list,
            comparison, comp_bambu, comp_FLAIR, comp_Lyric, comp_IsoTools, comp_Mandalorion,
            comp_Iso_IB, comp_FLAMES, comp_IsoQuant, comp_Spectra, comp_TALON_LAPA,
            comp_StringTie2, transcriptome_path2, annotation_path2, reference_path2,
            coverage2, coverage_dir2, cage_path2, quant_path2, poly_A_path2, platform2,
            library_preparation2, data_category2, sirv_list2, ercc_list2, sequin_list2
        ))
        thread.daemon = True
        thread.start()
        
        return jsonify({
            'success': True,
            'message': 'Your analysis has been started. You will receive an email when it completes.',
            'job_id': job.id
        })
        
    except Exception as e:
        print(f"Error starting job: {str(e)}", file=sys.stderr)
        return jsonify({'error': str(e)}), 500

def run_script_process_challenge1(job_id, job_results_dir, organism, transcriptome_path, annotation_path, reference_path,
                                 coverage, coverage_dir, cage_path, quant_path, poly_A_path, platform,
                                 library_preparation, tool, data_category, sirv_list, ercc_list, sequin_list,
                                 comparison, comp_bambu, comp_FLAIR, comp_Lyric, comp_IsoTools, comp_Mandalorion,
                                 comp_Iso_IB, comp_FLAMES, comp_IsoQuant, comp_Spectra, comp_TALON_LAPA,
                                 comp_StringTie2, transcriptome_path2, annotation_path2, reference_path2,
                                 coverage2, coverage_dir2, cage_path2, quant_path2, poly_A_path2, platform2,
                                 library_preparation2, data_category2, sirv_list2, ercc_list2, sequin_list2):
    """Background process for Challenge 1"""
    with app.app_context():
        job = Job.query.get(job_id)
        if not job:
            return
        
        logs_dir = "logs"
        os.makedirs(logs_dir, exist_ok=True)
        log_file_path = os.path.join(logs_dir, f"{job_id}.log")
        
        try:
            with open(log_file_path, 'w') as log_file:
                log_file.write(f"=== Challenge 1 Analysis Started ===\n")
                log_file.write(f"Job ID: {job_id}\n")
                log_file.write(f"Results Directory: {job_results_dir}\n")
                log_file.write(f"Started at: {datetime.utcnow().isoformat()}\n")
                log_file.write(f"Organism: {organism}\n")
                log_file.write(f"Platform: {platform}\n")
                log_file.write(f"Tool: {tool}\n")
                log_file.write(f"Library Preparation: {library_preparation}\n")
                log_file.write(f"Data Category: {data_category}\n")
                log_file.write(f"{'='*50}\n\n")
                log_file.flush()
                
                job.status = 'running'
                job.started_at = datetime.utcnow()
                db.session.commit()
                
                script_path = "lrgasp_event2_metrics/sqanti3_lrgasp.challenge1.py"
                
                # Build command arguments
                cmd_args = [
                    'python', script_path,
                    organism, transcriptome_path, annotation_path, reference_path, coverage, coverage_dir,
                    cage_path, quant_path, poly_A_path, platform, library_preparation, tool, data_category,
                    sirv_list, ercc_list, sequin_list, comparison, comp_bambu, comp_FLAIR, comp_Lyric,
                    comp_IsoTools, comp_Mandalorion, comp_Iso_IB, comp_FLAMES, comp_IsoQuant, comp_Spectra,
                    comp_TALON_LAPA, comp_StringTie2, transcriptome_path2, annotation_path2, reference_path2,
                    coverage2, coverage_dir2, cage_path2, quant_path2, poly_A_path2, platform2,
                    library_preparation2, data_category2, sirv_list2, ercc_list2, sequin_list2, '--gtf'
                ]
                
                # Add --dataset2 flag if second dataset is provided
                if transcriptome_path2 != 'NA':
                    cmd_args.insert(-1, '--dataset2')  # Insert before '--gtf'
                
                log_file.write(f"Executing command:\n{' '.join(cmd_args)}\n\n")
                log_file.write(f"{'='*50}\n")
                log_file.write(f"Script Output:\n")
                log_file.write(f"{'='*50}\n\n")
                log_file.flush()
                
                process = subprocess.Popen(
                    cmd_args,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True
                )
                
                for line in process.stdout:
                    line = line.rstrip()
                    print(line)
                    log_file.write(line + '\n')
                    log_file.flush()
                
                process.wait()
                
                log_file.write(f"\n{'='*50}\n")
                log_file.write(f"Script finished with return code: {process.returncode}\n")
                log_file.flush()
                
                if process.returncode != 0:
                    raise Exception(f"Script exited with code {process.returncode}")
                
                report_pattern = f"{job_results_dir}/results_file1/*.html"
                reports = glob.glob(report_pattern)
                if not reports:
                    raise Exception(f"No HTML report found in {job_results_dir}/results_file1/")
                
                source = reports[0]
                destination = f"static/report_{job_id}.html"
                shutil.copy(source, destination)
                
                log_file.write(f"Report copied from {source} to: {destination}\n")
                log_file.flush()
                
                job.status = 'completed'
                job.completed_at = datetime.utcnow()
                job.result_path = destination
                db.session.commit()
                
                log_file.write(f"\n{'='*50}\n")
                log_file.write(f"=== Challenge 1 Analysis Completed Successfully ===\n")
                log_file.write(f"Completed at: {datetime.utcnow().isoformat()}\n")
                log_file.write(f"{'='*50}\n")
                
                from utils import send_completion_email
                send_completion_email(job.id)
            
        except Exception as e:
            error_msg = f"Error in job {job_id}: {str(e)}"
            print(error_msg, file=sys.stderr)
            
            with open(log_file_path, 'a') as log_file:
                log_file.write(f"\n{'='*50}\n")
                log_file.write(f"=== ERROR ===\n")
                log_file.write(f"{error_msg}\n")
                log_file.write(f"{'='*50}\n")
            
            job.status = 'failed'
            job.completed_at = datetime.utcnow()
            job.error_message = str(e)
            db.session.commit()
            
            from utils import send_failure_email
            send_failure_email(job.id)

@app.route('/download_results/<int:job_id>')
def download_results(job_id):
    """Download results for a completed job"""
    job = Job.query.get_or_404(job_id)
    
    if job.status != 'completed' or not job.result_path:
        flash('Results are not available yet.', 'warning')
        return redirect(url_for('home'))
    
    if not os.path.exists(job.result_path):
        flash('Results file not found.', 'danger')
        return redirect(url_for('home'))
    
    return send_file(
        job.result_path,
        as_attachment=True,
        download_name=f'lrgasp_{job.challenge_type}_results_{job_id}.html'
    )

@app.route('/job_status/<int:job_id>')
def job_status(job_id):
    """Get status of a job"""
    job = Job.query.get_or_404(job_id)
    return jsonify({
        'id': job.id,
        'status': job.status,
        'challenge_type': job.challenge_type,
        'created_at': job.created_at.isoformat(),
        'completed_at': job.completed_at.isoformat() if job.completed_at else None,
        'error_message': job.error_message
    })

@app.route('/challenge_results')
def challenge_results():
    """Show challenge results page"""
    return render_template("show_report.html")

@app.route("/data")
def data():
    return render_template('data.html', title='Data')
