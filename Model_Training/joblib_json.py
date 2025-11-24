# fix_model_once_and_forget.py
import joblib
import xgboost as xgb

# Load your old joblib model
old_model = joblib.load(r"D:\MCED_DP\MCED_CANCER_DETECTOR.joblib")

# Save it in native XGBoost format (this format ALWAYS has predict_proba)
old_model.save_model(r"D:\MCED_DP\MCED_CANCER_DETECTOR.json")

print("DONE! Model saved as JSON → always works in Streamlit")