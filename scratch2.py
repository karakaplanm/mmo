import re

with open('/home/muka/mmo/week8.md', 'r') as f:
    lines = f.readlines()

new_lines = []
for i, line in enumerate(lines):
    if '```math' in line:
        # Check if previous line is not empty
        if len(new_lines) > 0 and new_lines[-1].strip() != '':
            new_lines.append('\n')
        new_lines.append(line)
    elif line.strip() == '```' and '```math' not in line:
        new_lines.append(line)
        # Check if next line is not empty (we just look ahead in lines)
        if i < len(lines) - 1 and lines[i+1].strip() != '':
            new_lines.append('\n')
    else:
        new_lines.append(line)

with open('/home/muka/mmo/week8.md', 'w') as f:
    f.writelines(new_lines)
