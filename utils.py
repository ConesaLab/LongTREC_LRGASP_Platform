import os
import sys
import subprocess
import shutil
from flask_mail import Mail, Message
from flask import current_app, url_for
from models import db, Job
from datetime import datetime
import zipfile

mail = Mail()

def send_email(to, subject, body, html=None):
    """Send email using Flask-Mail"""
    try:
        msg = Message(
            subject=subject,
            recipients=[to],
            body=body,
            html=html
        )
        mail.send(msg)
        return True
    except Exception as e:
        print(f"Error sending email: {str(e)}", file=sys.stderr)
        return False

def create_zip_file(source_dir, output_filename):
    """Create a zip file from a directory"""
    try:
        shutil.make_archive(output_filename.replace('.zip', ''), 'zip', source_dir)
        return True
    except Exception as e:
        print(f"Error creating zip file: {str(e)}", file=sys.stderr)
        return False

def send_completion_email(job_id):
    """Send completion email with download link"""
    with current_app.app_context():
        job = Job.query.get(job_id)
        if not job:
            print(f"Job {job_id} not found", file=sys.stderr)
            return
        
        # Create download URL
        download_url = f"{current_app.config['BASE_URL']}/download_results/{job.id}"
        
        subject = f"LRGASP {job.challenge_type.title()} - Analysis Complete"
        body = f"""
Hello,

Your LRGASP {job.challenge_type.title()} analysis has completed successfully!

You can download your results here:
{download_url}

This link will be available for 7 days.

Best regards,
LRGASP Platform Team
        """
        
        html = f"""
<html>
<body style="font-family: Arial, sans-serif; line-height: 1.6; color: #333;">
    <div style="max-width: 600px; margin: 0 auto; padding: 20px;">
        <h2 style="color: #198a82;">Analysis Complete!</h2>
        <p>Hello,</p>
        <p>Your LRGASP <strong>{job.challenge_type.title()}</strong> analysis has completed successfully!</p>
        <div style="margin: 30px 0; text-align: center;">
            <a href="{download_url}" 
               style="background-color: #198a82; color: white; padding: 12px 30px; 
                      text-decoration: none; border-radius: 5px; display: inline-block;">
                Download Results
            </a>
        </div>
        <p style="color: #666; font-size: 14px;">This link will be available for 7 days.</p>
        <hr style="border: none; border-top: 1px solid #ddd; margin: 30px 0;">
        <p style="color: #999; font-size: 12px;">
            Best regards,<br>
            LRGASP Platform Team
        </p>
    </div>
</body>
</html>
        """
        
        send_email(job.email, subject, body, html)

def send_failure_email(job_id):
    """Send failure email notification"""
    with current_app.app_context():
        job = Job.query.get(job_id)
        if not job:
            print(f"Job {job_id} not found", file=sys.stderr)
            return
        
        subject = f"LRGASP {job.challenge_type.title()} - Analysis Failed"
        body = f"""
Hello,

Unfortunately, your LRGASP {job.challenge_type.title()} analysis has failed.

Error: {job.error_message or 'Unknown error occurred'}

Please check your input files and try again. If the problem persists, contact support.

Best regards,
LRGASP Platform Team
        """
        
        html = f"""
<html>
<body style="font-family: Arial, sans-serif; line-height: 1.6; color: #333;">
    <div style="max-width: 600px; margin: 0 auto; padding: 20px;">
        <h2 style="color: #dc3545;">Analysis Failed</h2>
        <p>Hello,</p>
        <p>Unfortunately, your LRGASP <strong>{job.challenge_type.title()}</strong> analysis has failed.</p>
        <div style="background-color: #f8d7da; border: 1px solid #f5c6cb; border-radius: 5px; padding: 15px; margin: 20px 0;">
            <strong>Error:</strong> {job.error_message or 'Unknown error occurred'}
        </div>
        <p>Please check your input files and try again. If the problem persists, contact support.</p>
        <hr style="border: none; border-top: 1px solid #ddd; margin: 30px 0;">
        <p style="color: #999; font-size: 12px;">
            Best regards,<br>
            LRGASP Platform Team
        </p>
    </div>
</body>
</html>
        """
        
        send_email(job.email, subject, body, html)

def send_start_email(job_id):
    """Send email notification when analysis starts"""
    with current_app.app_context():
        job = Job.query.get(job_id)
        if not job:
            print(f"Job {job_id} not found", file=sys.stderr)
            return
        
        subject = f"LRGASP {job.challenge_type.title()} - Analysis Started"
        body = f"""
Hello,

Your LRGASP {job.challenge_type.title()} analysis has been started successfully!

Job ID: {job.id}
Started at: {job.created_at.strftime('%Y-%m-%d %H:%M:%S')}

You will receive another email when the analysis is complete with a link to download your results.

Best regards,
LRGASP Platform Team
        """
        
        html = f"""
<html>
<body style="font-family: Arial, sans-serif; line-height: 1.6; color: #333;">
    <div style="max-width: 600px; margin: 0 auto; padding: 20px;">
        <h2 style="color: #198a82;">Analysis Started!</h2>
        <p>Hello,</p>
        <p>Your LRGASP <strong>{job.challenge_type.title()}</strong> analysis has been started successfully!</p>
        <div style="background-color: #d1ecf1; border: 1px solid #bee5eb; border-radius: 5px; padding: 15px; margin: 20px 0;">
            <p style="margin: 5px 0;"><strong>Job ID:</strong> {job.id}</p>
            <p style="margin: 5px 0;"><strong>Started at:</strong> {job.created_at.strftime('%Y-%m-%d %H:%M:%S')}</p>
        </div>
        <p>You will receive another email when the analysis is complete with a link to download your results.</p>
        <hr style="border: none; border-top: 1px solid #ddd; margin: 30px 0;">
        <p style="color: #999; font-size: 12px;">
            Best regards,<br>
            LRGASP Platform Team
        </p>
    </div>
</body>
</html>
        """
        
        send_email(job.email, subject, body, html)

def allowed_file(filename, allowed_extensions):
    """Check if file has allowed extension"""
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in allowed_extensions

def evaluate_submission(file_path):
    # Placeholder for the actual evaluation logic
    # For example, read the file, compute metrics, and return the result
    # Here, we'll just return a dummy result
    return "Evaluation metrics would be displayed here."
