# Kinematics Results Database

import os
import sqlite3
import pandas as pd


class Database:

    def __init__(self, db_dir, db_name, table,
                 drop_on_init=False, columns=None):
        """
        kinematics propagation database class
        param db_dir: simulation results directory
        param db_name: database name i.e. filename
        param table: table name
        param columns: column names for table == output vars
        param drop_on_init: drop database if exists to avoid overlapping
        """
        self.table = table
        self.db = '%s%s.db' % (db_dir, db_name)
        if drop_on_init:
            self.drop()
        self.columns = columns
        self.conn = sqlite3.connect(self.db)
        self.c = self.conn.cursor()
        self.init_db()

    def init_db(self):
        """
        create simulation results database
        """
        self.c.execute("""CREATE TABLE IF NOT EXISTS {} ({})"""
                       .format(self.table, self.columns))
        print('...database ready\nusing %s table' % self.table)

    def insert(self, *args):
        """
        insert results into propagation database
        """
        n_cols = len(self.columns.split(','))
        n_vars = len(args[0])
        if n_cols != n_vars:
            raise Exception('DB and variable size mismatch')
            return
        mark = '?, ' * n_vars
        self.c.execute("""INSERT INTO {} ({}) VALUES ({})"""
                       .format(self.table, self.columns, mark[:-2]), (args[0]))
        self.conn.commit()

    def query(self, ti=None):
        """
        query table to retrieve results at time
        param ti: time to query
        If ti is None, read all results into pandas dataframe
        If ti is not None, return query results for that time
        return: pandas dataframe or row
        """
        if ti is None:
            cmd = "SELECT * FROM {}".format(self.table)
            self.c.execute(cmd)
            return pd.read_sql_query(cmd, self.conn)
        else:
            print('query for time, t = %f' % ti)
            cmd = "SELECT * FROM {} WHERE t={}".format(self.table, ti)
            self.c.execute(cmd)
            return self.c.fetchone()

    def drop(self):
        """
        drop filesystem database
        """
        try:
            os.remove(self.db)
            print('\n... %s database dropped\n' % self.db)
        except OSError:
            print('...database does not exist')

    def close(self):
        """
        close database connection
        """
        print('...closing database connection\n')
        self.conn.close()
