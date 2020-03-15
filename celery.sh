rabbitmq-server &
celery -A covid_profiler_web.worker worker --loglevel=info --concurrency=1

