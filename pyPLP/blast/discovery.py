import os
import shutil
import subprocess


def discover_blast_dbs(blastdb_dir: str) -> dict[str, str]:
    """Finds or builds available BLAST DB targets within the specified directory.

    Automatically handles:
    1. Subdirectories containing BLAST nucleotide indices (*.nin / *.nhr)
    2. Direct BLAST index files (*.nin / *.nhr) in blastdb_dir
    3. Unindexed organism FASTA files (*.fasta, *.fa, *.fna) in blastdb_dir or subdirectories,
       automatically running makeblastdb to index them on the fly.
    """
    if not os.path.isdir(blastdb_dir):
        return {}

    db_map = {}
    has_makeblastdb = shutil.which("makeblastdb") is not None

    for entry in sorted(os.listdir(blastdb_dir)):
        subpath = os.path.join(blastdb_dir, entry)

        # 1. Directory-based database
        if os.path.isdir(subpath):
            direct_prefix = os.path.join(subpath, entry)
            if os.path.exists(f"{direct_prefix}.nin") or os.path.exists(f"{direct_prefix}.nhr"):
                db_map[entry] = direct_prefix
                continue

            # Check for any *.nin / *.nhr in subpath
            found_nin = None
            for f in sorted(os.listdir(subpath)):
                if f.endswith((".nin", ".nhr")):
                    found_nin = os.path.join(subpath, os.path.splitext(f)[0])
                    break
            if found_nin:
                db_map[entry] = found_nin
                continue

            # Check for unindexed FASTA in subpath
            fasta_file = None
            for f in sorted(os.listdir(subpath)):
                if f.lower().endswith((".fasta", ".fa", ".fna")):
                    fasta_file = os.path.join(subpath, f)
                    break

            if fasta_file and has_makeblastdb:
                print(f"[BLAST DB Discovery] Indexing new organism '{entry}' in '{subpath}' using makeblastdb...")
                cmd = ["makeblastdb", "-in", fasta_file, "-dbtype", "nucl", "-out", direct_prefix, "-title", entry]
                res = subprocess.run(cmd, capture_output=True, text=True, check=False)
                if res.returncode == 0 and (os.path.exists(f"{direct_prefix}.nin") or os.path.exists(f"{direct_prefix}.nhr")):
                    print(f"  [+] Successfully indexed '{entry}' -> {direct_prefix}")
                    db_map[entry] = direct_prefix
                else:
                    print(f"[Warning] Failed to index '{fasta_file}': {res.stderr.strip() or res.stdout.strip()}")

        # 2. File placed directly in blastdb_dir -> auto-organize into its own subfolder
        elif os.path.isfile(subpath):
            name, ext = os.path.splitext(entry)
            ext = ext.lower()

            # Handle FASTA files or loose index files
            if ext in [".fasta", ".fa", ".fna", ".nin", ".nhr", ".ndb", ".njs", ".not", ".nsq", ".ntf", ".nto"]:
                org_dir = os.path.join(blastdb_dir, name)
                os.makedirs(org_dir, exist_ok=True)

                # Move all loose files sharing this organism stem into org_dir
                for loose in sorted(os.listdir(blastdb_dir)):
                    loose_path = os.path.join(blastdb_dir, loose)
                    if os.path.isfile(loose_path):
                        l_name, l_ext = os.path.splitext(loose)
                        if l_name == name:
                            dest = os.path.join(org_dir, loose)
                            if not os.path.exists(dest):
                                shutil.move(loose_path, dest)
                            elif loose_path != dest and os.path.exists(loose_path):
                                os.remove(loose_path)

                # Now check org_dir for existing index or build from fasta
                direct_prefix = os.path.join(org_dir, name)
                if os.path.exists(f"{direct_prefix}.nin") or os.path.exists(f"{direct_prefix}.nhr"):
                    db_map[name] = direct_prefix
                else:
                    fasta_file = None
                    for f in sorted(os.listdir(org_dir)):
                        if f.lower().endswith((".fasta", ".fa", ".fna")):
                            fasta_file = os.path.join(org_dir, f)
                            break
                    if fasta_file and has_makeblastdb:
                        print(f"[BLAST DB Discovery] Indexing new organism '{name}' in '{org_dir}' using makeblastdb...")
                        cmd = ["makeblastdb", "-in", fasta_file, "-dbtype", "nucl", "-out", direct_prefix, "-title", name]
                        res = subprocess.run(cmd, capture_output=True, text=True, check=False)
                        if res.returncode == 0 and (os.path.exists(f"{direct_prefix}.nin") or os.path.exists(f"{direct_prefix}.nhr")):
                            print(f"  [+] Successfully indexed '{name}' -> {direct_prefix}")
                            db_map[name] = direct_prefix
                        else:
                            print(f"[Warning] Failed to index '{fasta_file}': {res.stderr.strip() or res.stdout.strip()}")

    return db_map
