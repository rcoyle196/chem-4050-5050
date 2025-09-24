import pandas as pd

def write_dict_to_csv(x_values, data_dict, x_label="X", filename="output.csv"):
    data = {x_label: x_values}
    
    for key, values in data_dict.items():
        if len(values) != len(x_values):
            raise ValueError(f"Length mismatch for '{key}': expected {len(x_values)} values, got {len(values)}")
        data[key] = values
    
    df = pd.DataFrame(data)
    df.to_csv(filename, index=False)
    print(f"✅ Data saved to: {filename}")