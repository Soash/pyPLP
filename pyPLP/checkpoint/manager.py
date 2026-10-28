import datetime
import json
from pathlib import Path


class CheckpointManager:
    """Manages persistent step completion markers for pipeline resumption."""

    def __init__(self, job_dir: Path):
        self.job_dir = Path(job_dir)
        self.ckpt_dir = self.job_dir / ".checkpoints"
        self.ckpt_dir.mkdir(parents=True, exist_ok=True)

    def is_completed(self, step_name: str, expected_file: Path = None) -> bool:
        marker = self.ckpt_dir / f"{step_name}.done"
        if not marker.exists():
            return False
        if expected_file is not None and not Path(expected_file).exists():
            marker.unlink(missing_ok=True)
            return False
        return True

    def mark_completed(self, step_name: str, metadata: dict = None):
        marker = self.ckpt_dir / f"{step_name}.done"
        timestamp = datetime.datetime.now().isoformat()
        payload = {"step": step_name, "completed_at": timestamp}
        if metadata:
            payload.update(metadata)

        with open(marker, "w") as f:
            json.dump(payload, f, indent=2)
        print(f"[Checkpoint] Step '{step_name}' recorded -> {marker.name}")

    def get_metadata(self, step_name: str) -> dict:
        """Retrieves stored metadata for a completed checkpoint step."""
        marker = self.ckpt_dir / f"{step_name}.done"
        if marker.exists():
            try:
                with open(marker, "r") as f:
                    return json.load(f)
            except Exception:
                pass
        return {}

